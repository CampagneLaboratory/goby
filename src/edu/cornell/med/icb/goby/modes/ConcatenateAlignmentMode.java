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
import edu.cornell.med.icb.goby.aligners.AbstractAligner;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.alignments.processors.*;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.IOException;

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
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(ConcatenateAlignmentMode.class);

    public void setIgnoreTmh(boolean ignoreTmh) {
        this.ignoreTmh = ignoreTmh;
    }

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "concatenate-alignments";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Concatenate compact alignment files. Concatenation preserves sorting when " +
                    "every input alignment is already sorted.Reference sequences must match " +
                    "exactly across the input alignments.  Queries are assumed to be entirely " +
                    "distinct and will be treated as independent observations (e.g., reads from " +
                    "multiple independent samples). To this effect, alignment entries read from " +
                    "different input basenames, which would otherwise share an identical query " +
                    "index, are renumbered with distinct query indices (this behaviour can be " +
                    "turned off with the adjust-query-indices option).";

    private String[] inputFilenames;
    private String outputFile;
    /**
     * When this flag is true, the concatenation process does not produce a TMH file for the output alignment.
     */
    private boolean ignoreTmh;

    public void setAdjustSampleIndices(boolean adjustSampleIndices) {
        this.adjustSampleIndices = adjustSampleIndices;
    }

    private boolean adjustQueryIndices = true;
    private boolean realign = true;
    private AlignmentProcessorFactory alignmentProcessorFactory = new DefaultAlignmentProcessorFactory();
    private RandomAccessSequenceInterface genome;
    private boolean adjustSampleIndices;

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

    boolean upgrade = true;

    /**
     * Set this flag to false to prevent automatic alignment upgrading.
     *
     * @param upgrade
     */
    public void setUpgrade(boolean upgrade) {
        this.upgrade = upgrade;
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
        adjustSampleIndices = jsapResult.getBoolean("adjust-sample-indices", false);
        alignmentProcessorFactory = DiscoverSequenceVariantsMode.configureProcessor(jsapResult);
        genome = DiscoverSequenceVariantsMode.configureGenome(jsapResult);
        //final String codecName = jsapResult.getString("codec", null);

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
        final String[] basenames = AlignmentReaderImpl.getBasenames(inputFilenames);
        final boolean allSorted = isAllSorted(upgrade, basenames);
        if (allSorted) {
            System.out.println("input alignments are all sorted, the output will also be sorted.");
        } else {

            System.out.println("At least one of the input alignments is not sorted, the output will NOT be sorted.");

        }
        System.out.println(upgrade);
        // honor the no upgrade flag:
        final AlignmentReaderFactory alignmentReaderFactory = upgrade ?

                new DefaultAlignmentReaderFactory() :
                new NoUpgradeAlignmentReaderFactory();

        final ConcatAlignmentReader alignmentReader = allSorted ?
                new ConcatSortedAlignmentReader(alignmentReaderFactory, adjustQueryIndices, basenames) :
                new ConcatAlignmentReader(alignmentReaderFactory, adjustQueryIndices, basenames);

        alignmentReader.setAdjustSampleIndices(adjustSampleIndices);

        final ProgressLogger progress = new ProgressLogger();
        progress.displayFreeMemory = true;
        int entriesInOutputFile = 0;
        long numLogicalEntries = 0;
        long numEntries = 0;
        final int numQueries = alignmentReader.getNumberOfQueries();
        if (!ignoreTmh) {
            System.out.println("Concatenating TMH..");
            // too many hits is prepared as for Merge:
            Merge.prepareMergedTooManyHits(outputFile, alignmentReader.getNumberOfQueries(), 0, basenames);
            System.out.println("Done.");
        }
        progress.start("Concatenating entries");

        if (alignmentReader.getTargetLength() != null) {
            writer.setTargetLengths(alignmentReader.getTargetLength());
        }
        writer.setSorted(allSorted);


        AlignmentProcessorInterface processor = null;
        if (!allSorted) {
            processor = new DummyProcessorUnsorted(alignmentReader);
        } else {
            processor = alignmentProcessorFactory.create((ConcatSortedAlignmentReader) alignmentReader);
            if (processor instanceof LocalSortProcessor && genome == null) {
                System.err.println("A genome must be provided when realignment is requested.");
                System.exit(1);
            }
            processor.setGenome(genome);
        }
        assert processor != null : "processor cannot be null";
        Alignments.AlignmentEntry entry;
        while ((entry = processor.nextRealignedEntry(0, 0)) != null) {

            // query lengths are now always stored in the entry..
            writer.appendEntry(entry);

            numLogicalEntries += entry.getMultiplicity();
            numEntries++;
            entriesInOutputFile++;
            progress.lightUpdate();

        }
        alignmentReader.getStatistics();
        progress.done();

        writer.setSmallestSplitQueryIndex(alignmentReader.getSmallestSplitQueryIndex());
        writer.setLargestSplitQueryIndex(alignmentReader.getLargestSplitQueryIndex());

        if (alignmentReader.getQueryIdentifiers() != null) {
            writer.setQueryIdentifiers(alignmentReader.getQueryIdentifiers());
        }
        if (alignmentReader.getTargetIdentifiers() != null) {
            writer.setTargetIdentifiers(alignmentReader.getTargetIdentifiers());
        }
        writer.setNumQueries(alignmentReader.getNumberOfQueries());
        writer.setNumTargets(alignmentReader.getNumberOfTargets());

        writer.setStatistics(alignmentReader.getStatistics());
        writer.putStatistic("overall.matched.percent",
                String.format("%3.3g", divide(numLogicalEntries, numQueries) * 100d));
        writer.close();

        writer.printStats(System.out);
        System.out.printf("Wrote a total of %d alignment entries.%n", entriesInOutputFile);
        System.out.printf("Number of alignment entries realigned in the proximity of indels: %d (%3.3g %% of total)%n",
                processor.getModifiedCount(),
                divide(100 * processor.getModifiedCount(), processor.getProcessedCount()));

    }

    /**
     * Determine if all the alignments are sorted.
     *
     * @param basenames Basenames of a set of alignments.
     * @return True when all alignments are sorted, False otherwise.
     * @throws IOException If an error occurs reading an alignment.
     */
    public static boolean isAllSorted(final String[] basenames) throws IOException {
        return isAllSorted(true, basenames);
    }

    public static boolean isAllSorted(final boolean upgrade, final String[] basenames) throws IOException {
        boolean sorted = true;
        for (final String basename : basenames) {
            final AlignmentReader reader = new AlignmentReaderImpl(basename, upgrade);
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
