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
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.IterateAlignments;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import it.unimi.dsi.fastutil.longs.Long2LongMap;
import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongCollection;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

/**
 * Evaluate statistics for sequence variations found in alignments.
 *
 * @author Fabien Campagne
 */
public class SequenceVariationStatsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "sequence-variation-stats";

    /**
     * The overridden short mode name.
     */
    private static final String SHORT_MODE_NAME = "svs";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Evaluate statistics for sequence variations found in alignments";

    /**
     * The input filenames.
     */
    private String[] inputFilenames;

    /**
     * The output file.
     */
    private String outputFilename;
    /**
     * The input basenames.
     */
    private String[] basenames;
    private MyIterateAlignments alignmentIterator;


    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getShortModeName() {
        return SHORT_MODE_NAME;
    }

    enum OutputFormat {
        TSV,
        TAB_DELIMITED,

    }

    private OutputFormat outputFormat;

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

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
        basenames = AlignmentReaderImpl.getBasenames(inputFilenames);
        outputFilename = jsapResult.getString("output");
        outputFormat = OutputFormat.valueOf(jsapResult.getString("format").toUpperCase());

        alignmentIterator = new MyIterateAlignments();
        alignmentIterator.parseIncludeReferenceArgument(jsapResult);
        return this;
    }

    /**
     * Display sequence variations.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        PrintStream stream = null;
        try {
            stream = outputFilename == null ? System.out :
                    new PrintStream(new FileOutputStream(outputFilename));
            switch (outputFormat) {
                case TAB_DELIMITED:
                case TSV:
                    stream.println("basename\tread-index\tcount-variation-bases\tbases-at-index/all-variations-bases\tbases-at-index/all-reference-bases\tcount-reference-bases");
                    break;
            }

            for (final String basename : basenames) {
                alignmentIterator.resetTallies();
                final String[] singleBasename = {basename};
                // Iterate through each alignment and write sequence variations to output file:
                alignmentIterator.iterate(singleBasename);

                final Long2LongMap readIndexTallies = alignmentIterator.getReadIndexVariationTally();
                final double totalNumberOfVariationBases = sum(readIndexTallies.values());
                final double numberOfAlignmentEntries = alignmentIterator.getNumAlignmentEntries();

                final long countReferenceBases = alignmentIterator.getReferenceBaseCount();
                for (final long readIndex : readIndexTallies.keySet()) {
                    final long countVariationBasesAtReadIndex = readIndexTallies.get(readIndex);
                    final double frequency = ((double) countVariationBasesAtReadIndex) / totalNumberOfVariationBases;
                    final double alignFrequency = ((double) countVariationBasesAtReadIndex) / countReferenceBases;

                    stream.printf("%s\t%d\t%d\t%s\t%f\t%d%n",
                            FilenameUtils.getBaseName(basename),
                            readIndex,
                            countVariationBasesAtReadIndex,
                            frequency,
                            alignFrequency,
                            countReferenceBases);
                }
                stream.flush();
            }
        } finally {
            if (stream != System.out) {
                IOUtils.closeQuietly(stream);
            }
        }
    }

    private double sum(final LongCollection longCollection) {
        double sum = 0;
        for (final long value : longCollection) {
            sum += value;
        }
        return sum;
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
        new SequenceVariationStatsMode().configure(args).execute();
    }

    private static class MyIterateAlignments extends IterateAlignments {
        private Long2LongMap readIndexVariationTally = new Long2LongOpenHashMap();
        private final Long2LongMap readIndexReferenceTally = new Long2LongOpenHashMap();
        private long numAlignmentEntries;
        private long referenceBaseCount;

        public Long2LongMap getReadIndexVariationTally() {
            return readIndexVariationTally;
        }

        public long getNumAlignmentEntries() {
            return numAlignmentEntries;
        }

        @Override
        public void processAlignmentEntry(final AlignmentReader alignmentReader,
                                          final Alignments.AlignmentEntry alignmentEntry) {
            numAlignmentEntries += alignmentEntry.getMultiplicity();
            referenceBaseCount += alignmentEntry.getQueryLength();

            for (final Alignments.SequenceVariation variation : alignmentEntry.getSequenceVariationsList()) {

                final int toLength = variation.getTo().length();
                final int readIndexIncrementValue = alignmentEntry.getMatchingReverseStrand() ? -1 : 1;
                int readIndex = variation.getReadIndex();
                for (int i = 0; i < toLength; i++) {
                    if (readIndex < 0) {
                        System.err.printf("Negative read_index=%d for queryIndex=%d%n", readIndex, alignmentEntry.getQueryIndex());
                    }

                    final long value = readIndexVariationTally.get(readIndex);
                    final long changedBases = alignmentEntry.getMultiplicity();
                    readIndexVariationTally.put(readIndex, value + changedBases);

                    referenceBaseCount -= changedBases;
                    if (variation.getTo().charAt(i) != '-') {
                        readIndex += readIndexIncrementValue;
                    }
                }
            }
        }

        public void resetTallies() {
            readIndexVariationTally = new Long2LongOpenHashMap();
            numAlignmentEntries = 0;
            referenceBaseCount = 0;
        }

        public long getReferenceBaseCount() {
            return referenceBaseCount;
        }
    }
}
