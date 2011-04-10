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
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntCollection;
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

                final Int2IntMap readIndexTallies = alignmentIterator.getReadIndexVariationTally();
                final double totalNumberOfVariationBases = sum(readIndexTallies.values());
                final double numberOfAlignmentEntries = alignmentIterator.getNumAlignmentEntries();

                final long countReferenceBases = alignmentIterator.getReferenceBaseCount();
                for (final int readIndex : readIndexTallies.keySet()) {
                    final int countVariationBasesAtReadIndex = readIndexTallies.get(readIndex);
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

    private double sum(final IntCollection intCollection) {
        double sum = 0;
        for (final int value : intCollection) {
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
        private Int2IntMap readIndexVariationTally = new Int2IntOpenHashMap();
        private Int2IntMap readIndexReferenceTally = new Int2IntOpenHashMap();
        private int numAlignmentEntries;
        private int referenceBaseCount;

        public Int2IntMap getReadIndexVariationTally() {
            return readIndexVariationTally;
        }

        public int getNumAlignmentEntries() {
            return numAlignmentEntries;
        }

        @Override
        public void processAlignmentEntry(final AlignmentReaderImpl alignmentReader,
                                          final Alignments.AlignmentEntry alignmentEntry) {
            numAlignmentEntries += alignmentEntry.getMultiplicity();
            referenceBaseCount += alignmentEntry.getQueryLength();

            for (final Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {

                int toLength = var.getTo().length();
                int readIndex = var.getReadIndex();
                int readIndexIncrementValue = (alignmentEntry.getMatchingReverseStrand() ? -1 : 1);
                for (int i = 0; i < toLength; i++) {
                    if (readIndex < 0) {
                        System.err.printf("Negative read_index=%d for queryIndex=%d%n", readIndex, alignmentEntry.getQueryIndex());
                    }

                    final int value = readIndexVariationTally.get(readIndex);
                    final int changedBases = alignmentEntry.getMultiplicity();
                    readIndexVariationTally.put(readIndex, value + changedBases);

                    referenceBaseCount -= changedBases;
                    if (var.getTo().charAt(i) != '-') {
                        readIndex += readIndexIncrementValue;
                    }
                }
            }
        }

        public void resetTallies() {
            readIndexVariationTally = new Int2IntOpenHashMap();
            numAlignmentEntries = 0;
            referenceBaseCount = 0;
        }

        public int getReferenceBaseCount() {
            return referenceBaseCount;
        }
    }
}
