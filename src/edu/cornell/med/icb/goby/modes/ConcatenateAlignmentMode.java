/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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
import edu.cornell.med.icb.goby.aligners.AbstractAligner;
import edu.cornell.med.icb.goby.alignments.*;
import it.unimi.dsi.logging.ProgressLogger;

import java.io.File;
import java.io.IOException;
import java.io.FileNotFoundException;

/**
 * Concatenate compact alignment files. Concatenation preserves
 * sorting when every input alignment is sorted.
 * <p/>
 * Reference sequences must match exactly across the input alignments.
 * Query are assumed to be entirely distinct and will be treated as independent observations (e.g.,
 * reads from multiple independent samples). To this effect, alignment entries read from
 * different input basenames, which would otherwise share an identical query index,
 * are renumbered with distinct query indices.
 *
 * @author Fabien Campagne
 *         Date: Apr 28, 2009
 *         Time: 6:03:56 PM
 */
public class ConcatenateAlignmentMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "concatenate-alignments";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Concatenate compact alignment files. Concatenation preserves "
            + "sorting when every input alignment is sorted."
            + "Reference sequences must match exactly across the input alignments.  Queries "
            + "are assumed to be entirely distinct and will be treated as independent observations "
            + "(e.g., reads from multiple independent samples). To this effect, alignment entries "
            + "read from different input basenames, which would otherwise share an identical "
            + "query index, are renumbered with distinct query indices.";

    private String[] inputFilenames;
    private String outputFile;
    private boolean adjustQueryIndices = true;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    public void setAdjustQueryIndices(final boolean adjustQueryIndices) {
        this.adjustQueryIndices = adjustQueryIndices;
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

        outputFile = jsapResult.getString("output");
        adjustQueryIndices = jsapResult.getBoolean("adjust-query-indices", true);
        return this;
    }

    /**
     * Perform the concatenation.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {
        final String outputFilename = outputFile;
        final AlignmentWriter writer = new AlignmentWriter(outputFilename);
        final String[] basenames = AlignmentReader.getBasenames(inputFilenames);
        boolean allSorted = isAllSorted(basenames);
        if (allSorted) {
            System.out.println("input alignments are all sorted, the output will also be sorted.");
        } else {

            System.out.println("At least one of the input alignments is not sorted, the output will NOT be sorted.");

        }
        final ConcatAlignmentReader alignmentReader = allSorted ? new ConcatSortedAlignmentReader(adjustQueryIndices, basenames) :
                new ConcatAlignmentReader(adjustQueryIndices, basenames);
        final ProgressLogger progress = new ProgressLogger();
        progress.displayFreeMemory = true;
        int entriesInOutputFile = 0;
        long numLogicalEntries = 0;
        long numEntries = 0;
        final int numQueries = alignmentReader.getNumberOfQueries();

        progress.start("Concatenating entries");
        int[] queryLengths = null;
        if (alignmentReader.getQueryLengths() != null) {
            queryLengths = alignmentReader.getQueryLengths();
        }
        writer.setSorted(allSorted);
        for (final Alignments.AlignmentEntry entry : alignmentReader) {
            if (queryLengths != null) {
                writer.appendEntryWithLength(entry, queryLengths[entry.getQueryIndex()]);
            } else {
                writer.appendEntry(entry);
            }

            numLogicalEntries += entry.getMultiplicity();
            numEntries++;
            entriesInOutputFile++;
            progress.lightUpdate();

        }
        alignmentReader.getStatistics();
        progress.stop();
        // too many hits is prepared as for Merge:
        Merge.prepareMergedTooManyHits(outputFile, alignmentReader.getNumberOfQueries(), 0, basenames);

        writer.setNumTargets(alignmentReader.getNumberOfTargets());
        if (alignmentReader.getTargetIdentifiers() != null) {
            writer.setTargetIdentifiers(alignmentReader.getTargetIdentifiers());
        }
        if (alignmentReader.getTargetLength() != null) {
            writer.setTargetLengths(alignmentReader.getTargetLength());
        }
        if (alignmentReader.getQueryIdentifiers() != null) {
            writer.setQueryIdentifiers(alignmentReader.getQueryIdentifiers());
        }

        writer.setStatistics(alignmentReader.getStatistics());
        writer.putStatistic("overall.matched.percent",
                String.format("%3.3g", divide(numLogicalEntries, numQueries) * 100d));
        writer.close();

        writer.printStats(System.out);
        System.out.printf("Wrote a total of %d alignment entries.%n", entriesInOutputFile);
    }

    private boolean isAllSorted(String[] basenames) throws IOException {
        boolean sorted = true;
        for (String basename : basenames) {
            AlignmentReader reader = new AlignmentReader(basename);
            reader.readHeader();


            sorted &= reader.isSorted();

        }
        return sorted;

    }

    private double divide(final long a, final long b) {
        return ((double) a) / ((double) b);
    }


    public static void main(final String[] args) throws IOException, JSAPException {
        new ConcatenateAlignmentMode().configure(args).execute();
    }


    public void setInputFileNames(final String[] inputFileNames) {
        this.inputFilenames = inputFileNames;
    }

    public void setOutputFilename(final String outputFilename) {
        this.outputFile = outputFilename;
    }

    public File[] getOutputFiles() {
        return AbstractAligner.buildResults(outputFile);
    }
}
