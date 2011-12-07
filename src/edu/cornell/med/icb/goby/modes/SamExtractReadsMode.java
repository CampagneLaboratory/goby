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

import javax.management.RuntimeErrorException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Queue;
import java.util.concurrent.ArrayBlockingQueue;

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
     * The input filename.
     */
    private String inputFilename;

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

        inputFilename = jsapResult.getString("input");
        outputFilename = jsapResult.getString("output");

        if (!outputFilename.endsWith(".compact-reads")) {
            outputFilename = outputFilename + ".compact-reads";
        }
        qualityEncoding = QualityEncoding.valueOf(jsapResult.getString("quality-encoding").toUpperCase());
        processPairs = jsapResult.getBoolean("paired-end");
        processAtMost = jsapResult.getInt("process-at-most",Integer.MAX_VALUE);
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
            final ProgressLogger progress = new ProgressLogger(LOG);
            // the following is required to set validation to SILENT before loading the header (done in the SAMFileReader constructor)
            SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            final SAMFileReader parser = new SAMFileReader(new File(inputFilename), null);

            progress.start();
            int numReads = 0;
            Queue<SAMRecord> queue = new ArrayBlockingQueue<SAMRecord>(3);
            SAMRecordIterator samRecordIterator = parser.iterator();
            boolean hasMoreElements = samRecordIterator.hasNext();

            while (hasMoreElements) {
                if (queue.size() < 2) {
                    // enqueue
                    if (samRecordIterator.hasNext()) {
                        queue.add(samRecordIterator.next());
                    } else {

                        //we  may still have one read in the queue left to process:
                        hasMoreElements = !queue.isEmpty();
                    }
                }
                if (hasMoreElements && queue.size()>=2) {


                    // at least one element left, dequeue it
                    final SAMRecord first = queue.remove();

                    if (!processPairs) {

                        numReads = processSingleEndRead(writer, progress, numReads, first);
                    } else {
                        if (!queue.isEmpty()) {
                            final SAMRecord second = queue.remove();
                            if (first.getReadName().equals(second.getReadName())) {
                                numReads = processPairedEndRead(writer, progress, numReads, first, second);
                            } else {
                                System.err.printf("Error: when processing paired-end, two reads must have the same name in sequence. After #reads: %d Offending names %s %s %n",
                                        numReads,
                                        first.getReadName(), second.getReadName());
                                throw new RuntimeException();
                            }
                        } else {
                            System.err.printf("Error: when processing paired-end, two reads must be found for each name. After #reads: %d %n",
                                    numReads
                            );
                             throw new RuntimeException();
                        }
                    }
                }

                if (numReads > processAtMost) {
                    System.err.printf("Early stop after %d reads.%n", processAtMost);
                    break;
                }
            }
        } finally {
            writer.close();
        }
    }

    private int processPairedEndRead(ReadsWriter writer, ProgressLogger progress, int numReads, SAMRecord first,
                                     SAMRecord second) throws IOException {
        final String readId = first.getReadName();
        writer.setIdentifier(readId);
        writer.setSequence(byteToString(first.getReadBases()));

        writer.setQualityScores(FastaToCompactMode.convertQualityScores(qualityEncoding,
                byteToString(first.getBaseQualities()),
                false));
        writer.setPairSequence(byteToString(second.getReadBases()));

        writer.setQualityScoresPair(FastaToCompactMode.convertQualityScores(qualityEncoding,
                byteToString(second.getBaseQualities()),
                false));
        writer.appendEntry();
        progress.lightUpdate();
        numReads++;
        return numReads;
    }

    private int processSingleEndRead(ReadsWriter writer, ProgressLogger progress, int numReads, SAMRecord samRecord) throws IOException {
        final String readId = samRecord.getReadName();
        writer.setIdentifier(readId);
        writer.setSequence(byteToString(samRecord.getReadBases()));

        writer.setQualityScores(FastaToCompactMode.convertQualityScores(qualityEncoding,
                byteToString(samRecord.getBaseQualities()),
                false));
        writer.appendEntry();
        progress.lightUpdate();
        numReads++;
        return numReads;
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
