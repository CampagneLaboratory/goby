/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.reads.QualityEncoding;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Extract reads from a <a href="http://samtools.sourceforge.net/">SAM</a> file.
 *
 * @author Fabien Campagne
 */
public class SamExtractReadsMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(SamExtractReadsMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "sam-extract-reads";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Extract reads from a SAM file.";

    /**
     * The input filenames.
     */
    private String inputFilenames[];

    /**
     * The output file.
     */
    private String outputFilename;


    @Override
    public String getModeName() {
        return MODE_NAME;
    }


    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    private boolean processPairs = false;
    private String pairIndicator1 = null;
    private String pairIndicator2 = null;
    int processAtMost = -1;

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilenames = jsapResult.getStringArray("input");
        outputFilename = jsapResult.getString("output");
        if (!outputFilename.endsWith(".compact-reads")) {
            outputFilename = outputFilename + ".compact-reads";
        }
        qualityEncoding = QualityEncoding.valueOf(jsapResult.getString("quality-encoding").toUpperCase());

        processPairs = jsapResult.getBoolean("paired-end");
        final String tokens = jsapResult.getString("pair-indicator");
        if (processPairs && tokens != null) {
            final String[] tmp = tokens.split("[,]");
            if (tmp.length != 2) {
                System.err.println("Pair indicator argument must have exactly two tokens, separated by " +
                        "comma. Offending syntax: " + tokens);
                System.exit(1);
            }
            pairIndicator1 = tmp[0];
            pairIndicator2 = tmp[1];
        }
        processAtMost = jsapResult.getInt("process-at-most");
        return this;
    }

    private QualityEncoding qualityEncoding;

    /**
     * Display sequence variations.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        final ReadsWriter writer = new ReadsWriter(new FileOutputStream(outputFilename));
        try {
            for (String inputFilename : inputFilenames) {
                if (processPairs && inputFilename.contains(pairIndicator1)) {
                    processOneFile(inputFilename, writer);
                }
            }
        } finally {
            writer.close();
        }
    }

    private void processOneFile(String inputFilename, ReadsWriter writer) throws IOException {
        final ProgressLogger progress = new ProgressLogger(LOG);
        // the following is required to set validation to SILENT before loading the header (done in the SAMFileReader constructor)
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SAMFileReader parser = new SAMFileReader(new File(inputFilename), null);
        final String pairInputFilename;
        final SAMFileReader pairedParser;
        final SAMRecordIterator pairedIterator;
        progress.start();
        if (processPairs) {
            pairInputFilename = inputFilename.replace(pairIndicator1, pairIndicator2);
            LOG.info(String.format("Located paired-end input files (%s,%s)", inputFilename, pairInputFilename));
            pairedParser = new SAMFileReader(new File(pairInputFilename));
            pairedIterator = pairedParser.iterator();
        } else {
            pairedParser = null;
            pairedIterator = null;
            pairInputFilename = null;
        }
        int readNumber = 1;
        for (final SAMRecord samRecord : parser) {
            final String readId = samRecord.getReadName();
            writer.setIdentifier(readId);
            writer.setSequence(byteToString(samRecord.getReadBases()));

            writer.setQualityScores(FastaToCompactMode.convertQualityScores(qualityEncoding,
                    byteToString(samRecord.getBaseQualities()),
                    false));

            if (processPairs) {
                if (pairedIterator != null && !pairedIterator.hasNext()) {
                    System.err.printf("Paired file %s must have at least as many sequences as the primary input file %s. Aborting. %n",
                            inputFilename, pairInputFilename);
                    System.exit(1);
                } else {
                    assert pairedIterator != null;
                    final SAMRecord pairedSamRecord = pairedIterator.next();
                    if (!pairedSamRecord.getReadName().equals(readId)) {

                        System.err.printf("read identifier must match between paired file %s and primary input file %s." +
                                " Detected mismatch for read ids =%s %s at sequence number=%d. Aborting. %n",
                                inputFilename, pairInputFilename, readId, pairedSamRecord.getReadName(), readNumber);
                        System.exit(1);
                    }

                    writer.setPairSequence(byteToString(pairedSamRecord.getReadBases()));
                    writer.setQualityScoresPair(FastaToCompactMode.convertQualityScores(qualityEncoding,
                            byteToString(samRecord.getBaseQualities()),
                            false));
                }

            }
            readNumber++;
            if (processAtMost != -1 && readNumber > processAtMost) {
                break;
            }
            writer.appendEntry();
            progress.lightUpdate();
            parser.close();
            if (processPairs && pairedParser != null) {
                pairedParser.close();
            }
        }
    }


    private CharSequence byteToString(final byte[] input) {
        final MutableString buffer = new MutableString();
        buffer.setLength(input.length);
        for (int i = 0; i < input.length; i++) {
            buffer.setCharAt(i, (char) input[i]);

        }
        return buffer;
    }

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new SamExtractReadsMode().configure(args).execute();
    }
}
