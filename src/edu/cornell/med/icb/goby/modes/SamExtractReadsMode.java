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
import edu.cornell.med.icb.goby.alignments.EntryFlagHelper;
import edu.cornell.med.icb.goby.compression.MessageChunksWriter;
import edu.cornell.med.icb.goby.readers.sam.SAMRecordIterable;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.reads.QualityEncoding;
import edu.cornell.med.icb.goby.reads.ReadsWriterImpl;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionRegistry;
import it.unimi.dsi.Util;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.LinkedList;

/**
 * Extract reads from a <a href="http://samtools.sourceforge.net/">SAM</a> file.
 * WARNING: If the file contains pairs, the source SAM/BAM file >>MUST<< be sorted by read name.
 *
 * @author Fabien Campagne
 */
public class SamExtractReadsMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(SamExtractReadsMode.class);

    /**
     * Flag to indicate if log4j was configured.
     */
    private boolean debug;

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "sam-extract-reads";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Extract reads from a SAM file. WARNING: If the file contains pairs, the source SAM/BAM file >>MUST<< be sorted by read name.";

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
        processAtMost = jsapResult.getInt("process-at-most", Integer.MAX_VALUE);
        DynamicOptionRegistry.register(MessageChunksWriter.doc());

        return this;
    }

    private QualityEncoding qualityEncoding;

    private final LinkedList<SAMRecord> currentPairedReads = new LinkedList<SAMRecord>();
    private String currentPairedReadsName = null;

    private long numReadsTotal;
    private long numPairedReads;

    /**
     * Export reads from bam/sam to compact-reads
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        debug = Util.log4JIsConfigured();
        final ReadsWriter writer = new ReadsWriterImpl(new FileOutputStream(outputFilename));
        boolean finishEarly = false;

        numReadsTotal = 0;
        numPairedReads = 0;

        try {
            final ProgressLogger progress = new ProgressLogger(LOG);
            // the following is required to set validation to SILENT before loading the header (done in the SAMFileReader constructor)
            SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            final SAMFileReader parser = new SAMFileReader(new File(inputFilename), null);

            progress.start();

            for (final SAMRecord samRecord : new SAMRecordIterable(parser.iterator())) {
                if (!samRecord.getReadPairedFlag()) {
                    processSingleEndRead(writer, progress, samRecord);
                } else {
                    final String currentReadName = samRecord.getReadName();
                    if (currentPairedReadsName == null) {
                        // First paired read
                        currentPairedReadsName = currentReadName;
                        currentPairedReads.add(samRecord);
                    } else if (currentReadName.equals(currentPairedReadsName)) {
                        // Same read name as before
                        currentPairedReads.add(samRecord);
                    } else {
                        // New read name
                        processReadPairs(writer, progress, currentPairedReads);
                        currentPairedReads.add(samRecord);
                        currentPairedReadsName = currentReadName;
                    }
                }
                if (numReadsTotal > processAtMost) {
                    if (debug) {
                        LOG.info(String.format("Early stop after %d reads.", processAtMost));
                    }
                    finishEarly = true;
                    break;
                }
            }
            if (!finishEarly && !currentPairedReads.isEmpty()) {
                // End of file, output the file set of pairs if we have a set
                processReadPairs(writer, progress, currentPairedReads);
            }
            System.out.printf("Total number of reads written %d%n", numReadsTotal);
            System.out.printf("Number of paired-end reads written %d%n", numPairedReads);
        } finally {
            writer.close();
        }
    }

    private void processReadPairs(final ReadsWriter writer, final ProgressLogger progress,
                                 final LinkedList<SAMRecord> reads) throws IOException {
        if (reads.isEmpty()) {
            return;
        }
        String debugOuptut = null;
        if (debug && LOG.isDebugEnabled()) {
            final StringBuilder sb = new StringBuilder("\n");
            for (final SAMRecord read : reads) {
                sb.append(String.format("%s:[%s]%d->[%s]%d%n", read.getReadName(),
                        read.getReferenceName(), read.getAlignmentStart(),
                        read.getMateReferenceName(), read.getMateAlignmentStart()));
            }
            debugOuptut = sb.toString();
        }

        while (!reads.isEmpty()) {
            final SAMRecord samRecord1 = reads.removeFirst();
            final int mateAlignmentStart = samRecord1.getMateAlignmentStart();
            SAMRecord toRemove = null;
            for (final SAMRecord potentialMate : reads) {
                if (potentialMate.getAlignmentStart() == mateAlignmentStart) {
                    toRemove = potentialMate;
                    final SAMRecord first;
                    final SAMRecord second;
                    if (samRecord1.getFirstOfPairFlag()) {
                        first = samRecord1;
                        second = potentialMate;
                    } else {
                        first = potentialMate;
                        second = samRecord1;
                    }
                    processPairedEndRead(writer, progress, first, second);
                    break;
                }
            }
            if (toRemove != null) {
                reads.remove(toRemove);
                // We only need one successful pair pairing
                break;
            } else {
                if (reads.isEmpty()) {
                    processSingleEndRead(writer, progress, samRecord1);
                    if (debug && LOG.isDebugEnabled()) {
                        LOG.debug("Couldn't find a mate readName=" + samRecord1.getReadName() + ". Output will be non-paired. Flags=" + EntryFlagHelper.pairToString(samRecord1));
                        LOG.debug(debugOuptut);
                    }
                }
            }
        }

        reads.clear();
    }

    private void processPairedEndRead(final ReadsWriter writer, final ProgressLogger progress,
                                     final SAMRecord first, final SAMRecord second) throws IOException {
        final String readId = first.getReadName();
        writer.setIdentifier(readId);
        writer.setSequence(byteToString(buffer1, first.getReadBases()));

        writer.setQualityScores(FastaToCompactMode.convertQualityScores(qualityEncoding,
                byteToString(buffer2, first.getBaseQualities()),
                false));
        writer.setPairSequence(byteToString(buffer3, second.getReadBases()));

        writer.setQualityScoresPair(FastaToCompactMode.convertQualityScores(qualityEncoding,
                byteToString(buffer4, second.getBaseQualities()),
                false));
        writer.appendEntry();
        progress.lightUpdate();
        numPairedReads++;
        numReadsTotal++;
    }

    private void processSingleEndRead(final ReadsWriter writer, final ProgressLogger progress,
                                     final SAMRecord samRecord) throws IOException {
        final String readId = samRecord.getReadName();
        writer.setIdentifier(readId);
        writer.setSequence(byteToString(buffer1, samRecord.getReadBases()));

        writer.setQualityScores(FastaToCompactMode.convertQualityScores(qualityEncoding,
                byteToString(buffer2, samRecord.getBaseQualities()),
                false));
        writer.appendEntry();
        progress.lightUpdate();
        numReadsTotal++;
    }

    final MutableString buffer1 = new MutableString();
    final MutableString buffer2 = new MutableString();
    final MutableString buffer3 = new MutableString();
    final MutableString buffer4 = new MutableString();

    private CharSequence byteToString(final MutableString buffer, final byte[] input) {

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
        System.exit(0);
    }
}
