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
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math.stat.descriptive.rank.Percentile;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

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

    /** The minimum read length across all files. */
    private int minReadLength = Integer.MAX_VALUE;

    /** The maximum read length across all files. */
    private int maxReadLength = Integer.MIN_VALUE;

    /** The cumulative read length across all files. */
    private long cumulativeReadLength;

    /** The number of reads across all files. */
    private long numberOfReads;

    /** Whether or not to compute quantile information. */
    private boolean computeQuantiles;

    /** Number of quantiles used to characterize read length distribution. */
    private int numberOfQuantiles = 1;

    /** Display verbose output. */
    private boolean verbose;

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
     * @throws IOException error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        reset();
        final File[] inputFilesArray = jsapResult.getFileArray("input");
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
     * @throws IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        int numberOfFilesProcessed = 0;
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

        System.out.println();
        System.out.printf("Total number of files processed = %,d\n", numberOfFilesProcessed);
        System.out.printf("Total number of reads = %,d%n", numberOfReads);
        System.out.printf("Min read length = %,d%n", numberOfReads > 0 ? minReadLength : 0);
        System.out.printf("Max read length = %,d%n", numberOfReads > 0 ? maxReadLength : 0);
        System.out.printf("Avg read length = %,d%n", numberOfReads > 0 ? cumulativeReadLength / numberOfReads : 0);
        System.out.println();
    }

    /**
     * Print statistics about an alignment file in the Goby compact form.
     * @param file The file to display statistics about
     * @throws IOException if the file cannot be read
     */
    private void describeCompactAlignment(final File file) throws IOException {
        final String basename = AlignmentReader.getBasename(file.toString());
        System.out.printf("Compact Alignment basename = %s%n", basename);

        final AlignmentReader reader = new AlignmentReader(basename);
        reader.readHeader();
        System.out.println("Info from header:");
        System.out.printf("Number of query sequences = %,d%n", reader.getNumberOfQueries());
        System.out.printf("Number of target sequences = %,d%n", reader.getNumberOfTargets());
        System.out.printf("Has query identifiers = %s%n",
                reader.getQueryIdentifiers() != null && !reader.getTargetIdentifiers().isEmpty());
        System.out.printf("Has target identifiers = %s%n",
                reader.getTargetIdentifiers() != null && !reader.getTargetIdentifiers().isEmpty());

        int maxQueryIndex = -1;
        int maxTargetIndex = -1;
        int numEntries = 0;
        long numLogicalAlignmentEntries = 0;
        long total = 0;
        double avgScore = 0;
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
        }
        avgScore /= (double) numLogicalAlignmentEntries;

        System.out.printf("num query indices= %,d%n", maxQueryIndex + 1);
        System.out.printf("num target indices= %,d%n", maxTargetIndex + 1);
        System.out.printf("Number of alignment entries = %,d%n", numLogicalAlignmentEntries);
        System.out.printf("Percent matched = %3.2g%% %n",
                (double) numLogicalAlignmentEntries / (double) (long) maxQueryIndex * 100.0d);
        System.out.printf("Avg query alignment length = %,d%n", numEntries > 0 ? total / numEntries : -1);
        System.out.printf("Avg score alignment  = %f%n", avgScore);
    }

    /**
     * Print statistics about a reads file in the Goby compact form.
     * @param file The file to display statistics about
     * @throws IOException if the file cannot be read
     */
    private void describeCompactReads(final File file) throws IOException {
        System.out.printf("Compact reads filename = %s%n", file);

        // keep the read lengths for computing quantiles
        final DoubleArrayList readLengths = new DoubleArrayList();

        int minLength = Integer.MAX_VALUE;
        int maxLength = Integer.MIN_VALUE;

        int numberOfIdentifiers = 0;
        int numberOfDescriptions = 0;
        int numberOfSequences = 0;
        int numberOfQualityScores = 0;

        long numReadEntries = 0;
        long totalReadLength = 0;
        ReadsReader reader = null;
        try {
            reader = new ReadsReader(FileUtils.openInputStream(file));
            for (final Reads.ReadEntry entry : reader) {
                final int readLength = entry.getReadLength();

                // across this file
                numReadEntries++;
                totalReadLength += readLength;

                // across all files
                numberOfReads++;
                numberOfDescriptions += entry.hasDescription() ? 1 : 0;
                cumulativeReadLength += readLength;

                if (verbose && entry.hasDescription()) {
                    System.out.println("Description found: " + entry.getDescription());
                }
                numberOfIdentifiers += entry.hasReadIdentifier() ? 1 : 0;
                if (verbose && entry.hasReadIdentifier()) {
                    System.out.println("Identifier found: " + entry.getReadIdentifier());
                }
                numberOfSequences += entry.hasSequence() && !entry.getSequence().isEmpty() ? 1 : 0;
                numberOfQualityScores +=
                        entry.hasQualityScores() && !entry.getQualityScores().isEmpty() ? 1 : 0;

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
        } finally {
            if (reader != null) {
                reader.close();
            }
        }

        System.out.printf("Has identifiers    = %s (%,d) %n", numberOfIdentifiers > 0, numberOfIdentifiers);
        System.out.printf("Has descriptions   = %s (%,d) %n", numberOfDescriptions > 0, numberOfDescriptions);
        System.out.printf("Has sequences      = %s (%,d) %n", numberOfSequences > 0, numberOfSequences);
        System.out.printf("Has quality scores = %s (%,d) %n", numberOfQualityScores > 0, numberOfQualityScores);

        System.out.printf("Number of entries = %,d%n", numReadEntries);
        System.out.printf("Min read length = %,d%n", numReadEntries > 0 ? minLength : 0);
        System.out.printf("Max read length = %,d%n", numReadEntries > 0 ? maxLength : 0);
        System.out.printf("Avg read length = %,d%n", numReadEntries > 0 ? totalReadLength / numReadEntries : 0);

        // compute quantiles
        if (computeQuantiles) {
            final Percentile percentile = new Percentile();
            final double[] increasingReadLengths = readLengths.toDoubleArray();
            Arrays.sort(increasingReadLengths);
            System.out.printf("Read length quantiles = [ ");
            for (int quantile = 1; quantile < numberOfQuantiles + 1; quantile++) {
                System.out.printf("%,f ", percentile.evaluate(increasingReadLengths, quantile));
            }
            System.out.printf("]%n");
        }
    }

    /**
     * Get the maximum read length across all files processed so far.
     * @return The the maximum read length
     */
    public int getMaxReadLength() {
        return maxReadLength;
    }

    /**
     * Get the minimum read length across all files processed so far.
     * @return The the minimum read length
     */
    public int getMinReadLength() {
        return minReadLength;
    }

    /**
     * Get the cumulative length of the reads processed so far.
     * @return The total number of reads.
     */
    public long getCumulativeReadLength() {
        return cumulativeReadLength;
    }

    /**
     * Get the number of reads processed so far.
     * @return The total number of reads.
     */
    public long getNumberOfReads() {
        return numberOfReads;
    }

    /**
     * Get the list of files (reads/alignments) to process.
     * @return The list of files.
     */
    public List<File> getInputFiles() {
        return inputFiles;
    }

    /**
     * Add the specified file to the list of files to process.
     * @param inputFile The file to process
     */
    public void addInputFile(final File inputFile) {
        inputFiles.add(inputFile);
    }

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws JSAPException error parsing
     * @throws java.io.IOException error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new CompactFileStatsMode().configure(args).execute();
    }
}
