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
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentTooManyHitsReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.EntryFlagHelper;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import edu.cornell.med.icb.goby.algorithmic.data.DistinctIntValueCounterBitSet;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.apache.commons.math.stat.descriptive.rank.Percentile;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * Display some basic statistics on file in compact reads format.
 *
 * @author Kevin Dorff
 */
public class CompactFileStatsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "compact-file-stats";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Display some basic statistics on compact-reads and compact-alignment files.";

    /**
     * The input files.
     */
    private final List<File> inputFiles = new LinkedList<File>();

    /**
     * The minimum read length across all files.
     */
    private int minReadLength = Integer.MAX_VALUE;

    /**
     * The maximum read length across all files.
     */
    private int maxReadLength = Integer.MIN_VALUE;

    /**
     * The minimum quality length across all files.
     */
    private int minQualityLength = Integer.MAX_VALUE;

    /**
     * The maximum quality length across all files.
     */
    private int maxQualityLength = Integer.MIN_VALUE;

    /**
     * The cumulative read length across all files.
     */
    private long cumulativeReadLength;

    /**
     * The number of reads across all files.
     */
    private long numberOfReads;

    /**
     * Whether or not to compute quantile information.
     */
    private boolean computeQuantiles;

    /**
     * Number of quantiles used to characterize read length distribution.
     */
    private int numberOfQuantiles = 1;

    /**
     * Display verbose output.
     */
    private boolean verbose;

    /**
     * The output filename or null if we should write to stdout.
     */
    private String outputFilename;

    List<Boolean> pairedSamples = null;

    /**
     * The actual writer to use to write the output.
     */
    private PrintStream stream;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws IOException   error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        reset();
        final File[] inputFilesArray = jsapResult.getFileArray("input");
        outputFilename = jsapResult.getString("output");
        inputFiles.addAll(Arrays.asList(inputFilesArray));
        computeQuantiles = jsapResult.userSpecified("number-of-quantiles");
        numberOfQuantiles = jsapResult.getInt("number-of-quantiles", 1);
        verbose = jsapResult.getBoolean("verbose");
        return this;
    }

    /**
     * Reset the input file lists and cumulative statistics to default values.
     */
    public void reset() {
        inputFiles.clear();
        minReadLength = Integer.MAX_VALUE;
        maxReadLength = Integer.MIN_VALUE;
        numberOfReads = 0;
        cumulativeReadLength = 0;
    }

    /**
     * Run the FileStats mode.
     *
     * @throws IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        stream = null;
        try {
            stream = outputFilename == null ? System.out :
                    new PrintStream(new FileOutputStream(outputFilename));
            int numberOfFilesProcessed = 0;
            pairedSamples = new ArrayList<Boolean>(inputFiles.size());
            for (final File file : inputFiles) {
                if (file.exists() && file.canRead()) {
                    switch (FileExtensionHelper.determineCompactFileType(file)) {
                        case alignment:
                            describeCompactAlignment(file);
                            numberOfFilesProcessed++;
                            break;
                        case reads:
                            describeCompactReads(file);
                            numberOfFilesProcessed++;
                            break;
                        case unknown:
                        default:
                            System.err.println("Unknown file type: " + file);
                            break;
                    }
                } else {
                    System.err.println("Cannot read file: " + file);
                }
            }

            stream.println();
            stream.printf("Total number of files processed = %,d\n", numberOfFilesProcessed);
            stream.printf("All Compact-Reads files were paired-end = %s\n", isAllPairedSamples());
            stream.println();
        } finally {
            if (stream != System.out) {
                IOUtils.closeQuietly(stream);
            }
            stream = null;
        }
    }

    /**
     * Print statistics about an alignment file in the Goby compact form.
     *
     * @param file The file to display statistics about
     * @throws IOException if the file cannot be read
     */
    private void describeCompactAlignment(final File file) throws IOException {
        final String basename = AlignmentReader.getBasename(file.toString());
        stream.printf("Compact Alignment basename = %s%n", basename);

        final AlignmentReader reader = new AlignmentReader(basename);
        reader.readHeader();
        stream.println("Info from header:");
        stream.printf("Alignment produced by aligner=%s version=%s %n",reader.getAlignerName(),
                reader.getAlignerVersion());        
        stream.printf("Sorted: %b%n", reader.isSorted());
        stream.printf("Indexed: %b%n", reader.isIndexed());
        stream.printf("Number of target sequences = %,d%n", reader.getNumberOfTargets());
        final int[] targetLengthsFromHeader = reader.getTargetLength();
        stream.printf("Number of target length entries = %,d%n",
                ArrayUtils.getLength(reader.getTargetLength()));
        stream.printf("smallestSplitQueryIndex = %d%n",
                reader.getSmallestSplitQueryIndex());
        stream.printf("largestSplitQueryIndex = %d%n",
                reader.getLargestSplitQueryIndex());

        // simple statistics for target lengths
        final SummaryStatistics targetLengthStats = new SummaryStatistics();
        if (targetLengthsFromHeader != null) {
            for (final double d : targetLengthsFromHeader) {
                targetLengthStats.addValue(d);
            }
        }
        stream.printf("Min target length = %,d%n", (int) targetLengthStats.getMin());
        stream.printf("Max target length = %,d%n", (int) targetLengthStats.getMax());
        stream.printf("Mean target length = %,.2f%n", targetLengthStats.getMean());
        stream.println();

        stream.printf("Number of query sequences = %,d%n", reader.getNumberOfQueries());
        final int[] queryLengthsFromHeader = reader.getQueryLengths();
        stream.printf("Number of query length entries = %,d%n",
                ArrayUtils.getLength(queryLengthsFromHeader));

        final SummaryStatistics queryLengthStats = new SummaryStatistics();
        if (queryLengthsFromHeader != null) {
            stream.println("Query lengths are encoded in the header");
        }
        stream.println("Constant query lengths = " + reader.isConstantQueryLengths());

        stream.printf("Has query identifiers = %s%n",
                reader.getQueryIdentifiers() != null && !reader.getQueryIdentifiers().isEmpty());
        final IndexedIdentifier targetIdentifiers = reader.getTargetIdentifiers();
        final boolean hasTargetIdentifiers = targetIdentifiers != null && !targetIdentifiers.isEmpty();
        stream.printf("Has target identifiers = %s%n", hasTargetIdentifiers);
        if (verbose) {
            if (hasTargetIdentifiers) {
                for (Map.Entry<MutableString, Integer> entry : targetIdentifiers.entrySet()) {
                    stream.printf("  Target %s=%d with a length of %d%n", entry.getKey(), entry.getValue(), targetLengthsFromHeader[entry.getValue()]);
                }
            } else {
                for (Map.Entry<MutableString, Integer> entry : targetIdentifiers.entrySet()) {
                    stream.printf("  Target %d with a length of %d%n", entry.getValue(), targetLengthsFromHeader[entry.getValue()]);
                }
            }
        }
        stream.println();

        // the query indices that aligned. Includes those
        final DistinctIntValueCounterBitSet alignedQueryIndices = new DistinctIntValueCounterBitSet();

        describeAmbigousReads(basename, reader.getNumberOfQueries(), alignedQueryIndices);

        int maxQueryIndex = -1;
        int maxTargetIndex = -1;
        int numEntries = 0;
        long numLogicalAlignmentEntries = 0;
        long total = 0;
        double avgScore = 0;
        int sumNumVariations = 0;
        int numPaired = 0;
        int numProperlyPaired = 0;
        int numFirstInPair = 0;
        int numSecondInPair = 0;

        for (final Alignments.AlignmentEntry entry : reader) {
            numberOfReads++;   // Across all files
            numEntries++;      // Across this file
            numLogicalAlignmentEntries += entry.getMultiplicity();
            total += entry.getQueryAlignedLength();
            avgScore += entry.getScore();
            maxQueryIndex = Math.max(maxQueryIndex, entry.getQueryIndex());
            maxTargetIndex = Math.max(maxTargetIndex, entry.getTargetIndex());
            cumulativeReadLength += entry.getQueryAlignedLength();
            minReadLength = Math.min(minReadLength, entry.getQueryAlignedLength());
            maxReadLength = Math.max(maxReadLength, entry.getQueryAlignedLength());
            sumNumVariations += entry.getSequenceVariationsCount();
            alignedQueryIndices.observe(entry.getQueryIndex());

            // check entry then header for the query length
            if (entry.hasQueryLength()) {
                final double queryLength = entry.getQueryLength();
                queryLengthStats.addValue(queryLength);
            } else if (reader.hasQueryLengths() || reader.isConstantQueryLengths()) {
                final double queryLength = reader.getQueryLength(entry.getQueryIndex());
                queryLengthStats.addValue(queryLength);
            }
            numPaired += EntryFlagHelper.isPaired(entry) ? 1 : 0;
            numProperlyPaired += EntryFlagHelper.isProperlyPaired(entry) ? 1 : 0;
            numFirstInPair += EntryFlagHelper.isFirstInPair(entry) ? 1 : 0;
            numSecondInPair += EntryFlagHelper.isSecondInPair(entry) ? 1 : 0;
        }

        avgScore /= (double) numLogicalAlignmentEntries;

        final int numQuerySequences = reader.getNumberOfQueries();
        stream.printf("num query indices = %,d%n", numQuerySequences);
        final int numTargetSequences = maxTargetIndex + 1;
        final double avgNumVariationsPerQuery =
                ((double) sumNumVariations) / (double) numQuerySequences;
        stream.printf("num target indices = %,d%n", numTargetSequences);
        stream.printf("Number of alignment entries = %,d%n", numLogicalAlignmentEntries);
        stream.printf("Number of query indices that matched = %,d%n", alignedQueryIndices.count());
        stream.printf("Percent matched = %4.1f %% %n",
                (double) alignedQueryIndices.count() / (double) ((long) numQuerySequences) * 100.0d);
        stream.printf("Avg query alignment length = %,f%n",
                numEntries > 0 ? divide(total, numEntries) : -1);
        stream.printf("Avg score alignment = %f%n", avgScore);
        stream.printf("Avg number of variations per query sequence = %3.2f %n",
                avgNumVariationsPerQuery);
        final long size = file.length();
        stream.printf("Average bytes per entry = %f%n", divide(size, numLogicalAlignmentEntries));

        stream.printf("Min query length = %,d%n", (int) queryLengthStats.getMin());
        stream.printf("Max query length = %,d%n", (int) queryLengthStats.getMax());
        stream.printf("Mean query length = %,.2f%n", queryLengthStats.getMean());
        stream.printf("Percent paired reads = %,.2f %% %n", divide(numPaired, numQuerySequences * 2) * 100d);
        stream.printf("Percent properly paired reads = %,.2f %% %n", divide(numProperlyPaired, numQuerySequences * 2) * 100d);
        stream.printf("Percent first in pair = %,.2f %% %n", divide(numFirstInPair, numEntries) * 100d);
        stream.printf("Percent second in pair = %,.2f %% %n", divide(numSecondInPair, numEntries) * 100d);
    }

    private double divide(final long a, final long b) {
        return ((double) a) / (double) b;
    }

    private void describeAmbigousReads(final String basename, final double numReads, final DistinctIntValueCounterBitSet queryIndices) {
        try {
            final AlignmentTooManyHitsReader tmhReader = new AlignmentTooManyHitsReader(basename);
            queryIndices.observe(tmhReader.getQueryIndices());
            stream.printf("TMH: aligner threshold = %,d%n", tmhReader.getAlignerThreshold());
            stream.printf("TMH: number of ambiguous matches = %,d%n", tmhReader.getQueryIndices().size());
            stream.printf("TMH: %%ambiguous matches = %f %%%n", (tmhReader.getQueryIndices().size() * 100f) / numReads);
        } catch (IOException e) {
            stream.println("Cannot read TMH file for basename " + basename);
        }
    }

    /**
     * Print statistics about a reads file in the Goby compact form.
     *
     * @param file The file to display statistics about
     * @throws IOException if the file cannot be read
     */
    private void describeCompactReads(final File file) throws IOException {
        stream.printf("Compact reads filename = %s%n", file);

        // keep the read lengths for computing quantiles
        final DoubleArrayList readLengths = new DoubleArrayList();

        int minLength = Integer.MAX_VALUE;
        int maxLength = Integer.MIN_VALUE;

        int numberOfIdentifiers = 0;
        int numberOfDescriptions = 0;
        int numberOfSequences = 0;
        int numberOfSequencePairs = 0;
        int numberOfQualityScores = 0;
        int numberOfQualityScorePairs = 0;


        long totalReadLength = 0;
        long totalReadLengthPair = 0;
        final DistinctIntValueCounterBitSet allQueryIndices = new DistinctIntValueCounterBitSet();

        ReadsReader reader = null;
        boolean checkedForPaired = false;

        try {
            final long size = file.length();
            reader = new ReadsReader(FileUtils.openInputStream(file));
            for (final Reads.ReadEntry entry : reader) {
                final int readLength = entry.getReadLength();

                for (int i = 0; i < entry.getMetaDataCount(); i++) {
                    Reads.MetaData metaData = entry.getMetaData(i);
                    stream.printf("meta-data key=%s value=%s%n",
                            metaData.getKey(),
                            metaData.getValue());

                }

                // across this file
                allQueryIndices.observe(entry.getReadIndex());
                totalReadLength += readLength;
                totalReadLengthPair += entry.getReadLengthPair();

                // across all files
                numberOfReads++;
                numberOfDescriptions += entry.hasDescription() ? 1 : 0;
                cumulativeReadLength += readLength;

                if (verbose && entry.hasDescription()) {
                    stream.println("Description found: " + entry.getDescription());
                }
                numberOfIdentifiers += entry.hasReadIdentifier() ? 1 : 0;
                if (verbose && entry.hasReadIdentifier()) {
                    stream.println("Identifier found: " + entry.getReadIdentifier());
                }
                numberOfSequences += entry.hasSequence() && !entry.getSequence().isEmpty() ? 1 : 0;
                final boolean samplePaired = entry.hasSequencePair() && !entry.getSequencePair().isEmpty();
                if (samplePaired) {
                    numberOfSequencePairs += 1;
                }
                if (!checkedForPaired) {
                    // Check only the very first entry.
                    checkedForPaired = true;
                    pairedSamples.add(samplePaired);
                }
                if (entry.hasQualityScores() && !entry.getQualityScores().isEmpty()) {
                    numberOfQualityScores += 1;
                    final int qualityLength = entry.getQualityScores().size();
                    minQualityLength = Math.min(minQualityLength, qualityLength);
                    maxQualityLength = Math.max(maxQualityLength, qualityLength);
                }

                numberOfQualityScorePairs +=
                        entry.hasQualityScoresPair() && !entry.getQualityScoresPair().isEmpty() ? 1 : 0;

                // we only need to keep all the read lengths if quantiles are being computed
                if (computeQuantiles) {
                    readLengths.add(readLength);
                }
                minLength = Math.min(minLength, readLength);
                maxLength = Math.max(maxLength, readLength);

                // adjust the min/max length of across all files
                minReadLength = Math.min(minReadLength, readLength);
                maxReadLength = Math.max(maxReadLength, readLength);
            }


            stream.printf("Average bytes per entry: %f%n", divide(size, allQueryIndices.count()));
            stream.printf("Average bytes per base: %f%n", divide(size, cumulativeReadLength));
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
        final int numReadEntries = allQueryIndices.count();
        stream.printf("Has identifiers = %s (%,d) %n", numberOfIdentifiers > 0, numberOfIdentifiers);
        stream.printf("Has descriptions = %s (%,d) %n", numberOfDescriptions > 0, numberOfDescriptions);
        stream.printf("Has sequences = %s (%,d) %n", numberOfSequences > 0, numberOfSequences);
        stream.printf("Has sequencePairs = %s (%,d) %n", numberOfSequencePairs > 0, numberOfSequencePairs);
        stream.printf("Has quality scores = %s (%,d) %n", numberOfQualityScores > 0, numberOfQualityScores);
        stream.printf("Has quality score Pairs = %s (%,d) %n", numberOfQualityScorePairs > 0, numberOfQualityScorePairs);

        stream.printf("Number of entries = %,d%n", numReadEntries);
        stream.printf("Min read length = %,d%n", numReadEntries > 0 ? minLength : 0);
        stream.printf("Max read length = %,d%n", numReadEntries > 0 ? maxLength : 0);
        stream.printf("Min quality length = %,d%n", numberOfQualityScores > 0 ? minQualityLength : 0);
        stream.printf("Max quality length = %,d%n", numberOfQualityScores > 0 ? maxQualityLength : 0);
        stream.printf("Avg read length = %,d%n", numReadEntries > 0 ? totalReadLength / numReadEntries : 0);
        stream.printf("Avg read pair length = %,d%n", numReadEntries > 0 ? totalReadLengthPair / numReadEntries : 0);

        // compute quantiles
        if (computeQuantiles) {
            final Percentile percentile = new Percentile();
            final double[] increasingReadLengths = readLengths.toDoubleArray();
            Arrays.sort(increasingReadLengths);
            stream.printf("Read length quantiles = [ ");
            for (int quantile = 1; quantile < numberOfQuantiles + 1; quantile++) {
                stream.printf("%,f ", percentile.evaluate(increasingReadLengths, quantile));
            }
            stream.printf("]%n");
        }
    }

    /**
     * Get the maximum read length across all files processed so far.
     *
     * @return The the maximum read length
     */
    public int getMaxReadLength() {
        return maxReadLength;
    }

    /**
     * Get the minimum read length across all files processed so far.
     *
     * @return The the minimum read length
     */
    public int getMinReadLength() {
        return minReadLength;
    }

    /**
     * Get if all the files scanned contains paired samples (only checks the FIRST sample in each file).
     *
     * @return if all the files contained paired samples
     */
    public boolean isAllPairedSamples() {
        boolean allPairedSamples = true;
        if (pairedSamples == null || pairedSamples.size() == 0) {
            allPairedSamples = false;
        } else {
            for (boolean pairedSample : pairedSamples) {
                if (!pairedSample) {
                    allPairedSamples = false;
                    break;
                }
            }
        }
        return allPairedSamples;
    }

    /**
     * A list if which of the scanned samples were paired and which weren't. Note that
     * the pairing scan only checks the FIRST sample in the reads file.
     *
     * @return the list of which scanned samples were paired (one per file scanned).
     */
    public List<Boolean> getPairedSamples() {
        return pairedSamples;
    }

    /**
     * Get the cumulative length of the reads processed so far.
     *
     * @return The total number of reads.
     */
    public long getCumulativeReadLength() {
        return cumulativeReadLength;
    }

    /**
     * Get the number of reads processed so far.
     *
     * @return The total number of reads.
     */
    public long getNumberOfReads() {
        return numberOfReads;
    }

    /**
     * Get the list of files (reads/alignments) to process.
     *
     * @return The list of files.
     */
    public List<File> getInputFiles() {
        return inputFiles;
    }

    /**
     * Add the specified file to the list of files to process.
     *
     * @param inputFile The file to process
     */
    public void addInputFile(final File inputFile) {
        inputFiles.add(inputFile);
    }

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws JSAPException       error parsing
     * @throws java.io.IOException error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new CompactFileStatsMode().configure(args).execute();
    }
}
