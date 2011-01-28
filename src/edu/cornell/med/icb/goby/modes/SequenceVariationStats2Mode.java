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

import java.io.IOException;
import java.io.PrintStream;
import java.io.FileOutputStream;

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.IterateSortedAlignments;
import edu.cornell.med.icb.goby.alignments.ConcatSortedAlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;

/**
 * @author Fabien Campagne
 *         Date: Sep 8, 2010
 *         Time: 12:42:59 PM
 */
/**
 * Evaluate statistics for sequence variations found in alignments. Alternative implementation with the IterateSortedAlignmentsHelper.
 *
 * @author Fabien Campagne
 */
public class SequenceVariationStats2Mode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "sequence-variation-stats2";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Evaluate statistics for sequence variations found in alignments. (alternative implementation.)";

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
        basenames = AlignmentReader.getBasenames(inputFilenames);
        outputFilename = jsapResult.getString("output");
        outputFormat = OutputFormat.valueOf(jsapResult.getString("format").toUpperCase());


        this.jsapResult = jsapResult;
        return this;
    }

    JSAPResult jsapResult;

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
                    stream.println("basename\tread-index\tcount-variation-bases\tbases-at-index/all-variations-bases\tbases-at-index/all-reference-bases\tcount-reference-bases\tcount-reference-bases-at-index");
                    break;
            }


            for (final String basename : basenames) {
                MyIterateSortedAlignments iterator = new MyIterateSortedAlignments();
                iterator.parseIncludeReferenceArgument(jsapResult);

                final String[] singleBasename = {basename};
                // Iterate through each alignment and write sequence variations to output file:
                iterator.iterate(singleBasename);

                final int[] readIndexVariationTallies = iterator.getReadIndexVariationTally();
                final int[] readIndexReferenceTallies = iterator.getReadIndexReferenceTally();
                final double totalNumberOfVariationBases = sum(readIndexVariationTallies);
                final double numberOfAlignmentEntries = iterator.getNumAlignmentEntries();

                final long countReferenceBases = iterator.getReferenceBaseCount();
                int maxReadIndex = iterator.getMaxReadIndex();
                for (int readIndex = 1; readIndex < maxReadIndex + 1; readIndex++) {
                    final int countVariationBasesAtReadIndex = readIndexVariationTallies[readIndex];
                    final int countReferenceBasesAtReadIndex = readIndexReferenceTallies[readIndex];
                    final double frequency = ((double) countVariationBasesAtReadIndex) / totalNumberOfVariationBases;
                    final double alignFrequency = ((double) countVariationBasesAtReadIndex) / countReferenceBases;

                    stream.printf("%s\t%d\t%d\t%s\t%f\t%d\t%d%n",
                            FilenameUtils.getBaseName(basename),
                            readIndex,
                            countVariationBasesAtReadIndex,
                            frequency,
                            alignFrequency,
                            countReferenceBases,
                            countReferenceBasesAtReadIndex);
                }
                stream.flush();
            }
        } finally {
            if (stream != System.out) {
                IOUtils.closeQuietly(stream);
            }
        }
    }

    private double sum(final int[] intCollection) {
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


    private class CountsAtPosition {
    }

    private class MyIterateSortedAlignments extends IterateSortedAlignments<CountsAtPosition> {
        private int[] readIndexVariationTally;
        private int[] readIndexReferenceTally;
        private int referenceBaseCount;
        private int maxReadIndex = -1;

        private MyIterateSortedAlignments() {
            readIndexReferenceTally = new int[100000];
            readIndexVariationTally = new int[100000];
        }

        public void observeReferenceBase(ConcatSortedAlignmentReader sortedReaders,
                                         Alignments.AlignmentEntry alignmentEntry,
                                         Int2ObjectMap<CountsAtPosition> positionToBases,
                                         int queryLength, boolean forwardStrand,
                                         int positionInMatch) {

        }

        public void observeReferenceBase(ConcatSortedAlignmentReader sortedReaders, Alignments.AlignmentEntry alignmentEntry,
                                         Int2ObjectMap<CountsAtPosition> positionToBases,
                                         int currentReferenceIndex, int currentRefPosition, int readIndex) {
            if (readIndex >= 1) {
         //       assert readIndex > 0 : String.format("positionInMatch=%d %s %n", positionInMatch, alignmentEntry);
                maxReadIndex = Math.max(maxReadIndex, readIndex);
                int count = readIndexReferenceTally[readIndex];
                readIndexReferenceTally[readIndex] = count + 1;
                referenceBaseCount += 1;
            }
        }

        public void observeVariantBase(ConcatSortedAlignmentReader sortedReaders,
                                       Alignments.AlignmentEntry alignmentEntry, Int2ObjectMap<CountsAtPosition> positionToBases,
                                       Alignments.SequenceVariation var,
                                       char toChar, char fromChar, int currentReferenceIndex, int currentRefPosition, int currentReadIndex) {

            maxReadIndex = Math.max(maxReadIndex, currentReadIndex);
            int count = readIndexVariationTally[currentReadIndex];
            readIndexVariationTally[currentReadIndex] = count + 1;
        }

        public void processPositions(int referenceIndex, int intermediatePosition, CountsAtPosition positionBaseInfos) {

        }

        public int getReferenceBaseCount() {
            return referenceBaseCount;
        }

        public int[] getReadIndexVariationTally() {
            return readIndexVariationTally;
        }

        public int[] getReadIndexReferenceTally() {
            return readIndexReferenceTally;
        }


        public int getMaxReadIndex() {
            return maxReadIndex;
        }
    }
}
