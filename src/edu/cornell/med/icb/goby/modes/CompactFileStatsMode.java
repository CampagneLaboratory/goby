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
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import org.apache.commons.math.stat.descriptive.rank.Percentile;

import java.io.File;
import java.io.FileInputStream;
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
    public static final String MODE_NAME = "compact-file-stats";
    public static final String MODE_DESCRIPTION =
            "Display some basic statistics on compact-reads and compact-alignment files.";

    /**
     * The input files.
     */
    private final List<File> inputFiles = new LinkedList<File>();

    /** The min read length. */
    private int minReadLength = Integer.MAX_VALUE;

    /** The max read length. */
    private int maxReadLength = Integer.MIN_VALUE;

    private int numberOfQuantiles = 1;

    /** The number of reads. */
    private long numberOfReads;

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
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        reset();
        final File[] inputFilesArray = jsapResult.getFileArray("input");
        inputFiles.addAll(Arrays.asList(inputFilesArray));
        numberOfQuantiles = jsapResult.getInt("number-of-quantiles", 1);
        return this;
    }

    public void reset() {
        inputFiles.clear();
        minReadLength = Integer.MAX_VALUE;
        maxReadLength = Integer.MIN_VALUE;
        numberOfReads = 0;
    }

    /**
     * Run the FileStats mode.
     * @throws IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        for (final File file : inputFiles) {
            if (file.exists() && file.canRead()) {
                switch (FileExtensionHelper.determineCompactFileType(file)) {
                    case alignment:
                        describeCompactAlignment(file.toString());
                        break;
                    case reads:
                        describeCompactReads(file.toString());
                        break;
                    case unknown:
                        System.err.println("Unknown file type: " + file);
                        break;
                }
            } else {
                System.err.println("Cannot read file: " + file);
            }
        }
    }

    private void describeCompactAlignment(final String file) throws IOException {
        final String basename = AlignmentReader.getBasename(file);
        System.out.printf("Compact Alignment basename = %s%n", basename);

        final AlignmentReader reader = new AlignmentReader(basename);
        reader.readHeader();
        System.out.println("Info from header:");
        System.out.printf("Number of query sequences = %d%n", reader.getNumberOfQueries());
        System.out.printf("Number of target sequences = %d%n", reader.getNumberOfTargets());
        System.out.printf("has query identifiers = %s%n",
                reader.getQueryIdentifiers() != null && !reader.getTargetIdentifiers().isEmpty());
        System.out.printf("has target identifiers = %s%n",
                reader.getTargetIdentifiers() != null && !reader.getTargetIdentifiers().isEmpty());

        int maxQueryIndex = -1;
        int maxTargetIndex = -1;
        int numEntries = 0;
        long numLogicalAlignmentEntries = 0;
        long total = 0;
        double avgScore = 0;
        for (final Alignments.AlignmentEntry entry : reader) {
            numberOfReads++; // Across all files
            numEntries++; // Across this file
            numLogicalAlignmentEntries+=entry.getMultiplicity();
            total += entry.getQueryAlignedLength();
            avgScore += entry.getScore();
            maxQueryIndex = Math.max(maxQueryIndex, entry.getQueryIndex());
            maxTargetIndex = Math.max(maxTargetIndex, entry.getTargetIndex());
            minReadLength = Math.min(minReadLength, entry.getQueryAlignedLength());
            maxReadLength = Math.max(maxReadLength, entry.getQueryAlignedLength());
        }
        avgScore /= (double) numLogicalAlignmentEntries;

        System.out.printf("num query indices= %d%n", maxQueryIndex + 1);
        System.out.printf("num target indices= %d%n", maxTargetIndex + 1);
        System.out.printf("Number of alignment entries = %d%n", numLogicalAlignmentEntries);
        System.out.printf("Percent matched = %3.2g%% %n", divide(numLogicalAlignmentEntries, maxQueryIndex) * 100.0d);
        System.out.printf("Avg query alignment length = %,d%n", numEntries>0 ? total / numEntries : -1);
        System.out.printf("Avg score alignment  = %f%n", avgScore);
    }

    private double divide(final long a, final long b) {
        return  (double)a / (double)b;
    }

    private void describeCompactReads(final String file) throws IOException {
        final DoubleArrayList readLengths = new DoubleArrayList();
        System.out.printf("Compact reads filename = %s%n", file);

        int numberOfIdentifiers = 0;
        int numberOfDescriptions = 0;
        int numberOfSequences = 0;

        long numReadEntries = 0;
        long totalReadLength = 0;
        ReadsReader reader = null;
        try {
            reader = new ReadsReader(new FileInputStream(file));
            for (final Reads.ReadEntry entry : reader) {
                numberOfReads++;  // across all files
                numReadEntries++; // across this file
                totalReadLength += entry.getReadLength();
                numberOfDescriptions += entry.hasDescription() ? 1 : 0;
                if (entry.hasDescription()) {
                    System.out.println("Description found: " + entry.getDescription());
                }
                numberOfIdentifiers += entry.hasReadIdentifier() ? 1 : 0;
                if (entry.hasReadIdentifier()) {
                    System.out.println("Identifier found: " + entry.getReadIdentifier());
                }
                numberOfSequences += entry.hasSequence() && !entry.getSequence().isEmpty() ? 1 : 0;
                final int readLength = entry.getReadLength();
                readLengths.add(readLength);
                minReadLength = Math.min(minReadLength, readLength);
                maxReadLength = Math.max(maxReadLength, readLength);
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }

        System.out.printf("has identifiers  = %s (%d) %n", numberOfIdentifiers > 0, numberOfIdentifiers);
        System.out.printf("has descriptions = %s (%d) %n", numberOfDescriptions > 0, numberOfDescriptions);
        System.out.printf("has sequences    = %s (%d) %n", numberOfSequences > 0, numberOfSequences);

        System.out.printf("Number of entries = %d%n", numReadEntries);
        System.out.printf("Min read length = %,d%n", numReadEntries > 0 ? minReadLength : 0);
        System.out.printf("Max read length = %,d%n", numReadEntries > 0 ? maxReadLength : 0);
        System.out.printf("Avg read length = %,d%n", numReadEntries > 0 ? totalReadLength / numReadEntries : 0);

        // compute quantiles
        final Percentile percentile = new Percentile();
        final double[] increasingReadLengths = readLengths.toDoubleArray();
        Arrays.sort(increasingReadLengths);
        System.out.printf("Read length quantiles = [ ");
        for (int quantile = 1; quantile < numberOfQuantiles + 1; quantile++) {
            System.out.printf("%,f ", percentile.evaluate(increasingReadLengths, quantile));
        }
        System.out.printf("]%n");
    }

    public int getMaxReadLength() {
        return maxReadLength;
    }

    public int getMinReadLength() {
        return minReadLength;
    }

    public long getNumberOfReads() {
        return numberOfReads;
    }

    public List<File> getInputFiles() {
        return inputFiles;
    }

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
