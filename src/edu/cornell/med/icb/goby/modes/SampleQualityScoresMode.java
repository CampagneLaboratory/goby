/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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
import edu.cornell.med.icb.goby.readers.FastXEntry;
import edu.cornell.med.icb.goby.readers.FastXReader;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Read the first n (defaults to 10,000) entries of a compact-reads or FASTQ file and report the minimum,
 * maximum, and average quality score value and provide a guess at the quality encoding scheme
 * (Phred, Sanger, Illumina/Solexa). It's worth noting that compact-reads files, if created correctly,
 * should always report Phred.
 *
 * @author Kevin C. Dorff
 *         Date: Sept 13, 2011
 */
public class SampleQualityScoresMode extends AbstractGobyMode {
    /**
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(SampleQualityScoresMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "sample-quality-scores";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Read the first n (defaults to 10,000) entries of a compact-reads or FASTQ file and report the minimum, " +
            "maximum, and average quality score value and provide a guess at the quality encoding scheme " +
            "(Phred, Sanger, Illumina/Solexa). It's worth noting that compact-reads files, if created correctly, " +
            "should always report Phred. ";

    /**
     * The reads files to process.
     */
    private final List<String> inputFilenames = new ArrayList<String>();

    /**
     * The CURRENT min found qual score (only the last file processed).
     * See avgQualScores and likelyEncodings file results.
     */
    private int minQualScore;

    /**
     * The CURRENT max found qual score (only the last file processed).
     * See avgQualScores and likelyEncodings file results.
     */
    private int maxQualScore;

    /**
     * The CURRENT sum of all of the quality scores, so we can perform an average (only the last file processed).
     * See avgQualScores and likelyEncodings file results.
     */
    private long sumQualScores;

    /**
     * The CURRENT number of qual scores sampled, so we can perform an average (only the last file processed).
     * See avgQualScores and likelyEncodings file results.
     */
    private long numQualScoresSampled;

    /**
     * The average qual scores for each processed file.
     */
    private final IntList avgQualScores = new IntArrayList();

    /**
     * The likely encodings for each processed file.
     */
    private final List<String> likelyEncodings = new ArrayList<String>();

    /**
     * The number of read entries from each input file to to process.
     */
    private int numberOfReadEntriesToProcess = 10000;

    /**
     * If quality scores were found in the current input files.
     */
    boolean qualityScoresFound;
    
    /**
     * Returns the mode name defined by subclasses.
     *
     * @return The name of the mode
     */
    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    /**
     * Returns the mode description defined by subclasses.
     *
     * @return A description of the mode
     */
    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Clear the list of files to process.
     *
     * @param inputFilename The list of files.
     */
    public void addInputFilename(final String inputFilename) {
        inputFilenames.add(inputFilename);
    }

    /**
     * Clear the list of files to process.
     */
    public void clearInputFilenames() {
        inputFilenames.clear();
    }

    /**
     * The list of input filenames.
     * @return list of input filenames.
     */
    public List<String> getInputFilenames() {
        return inputFilenames;
    }

    /**
     * Get the number of read entries to process per input file.
     *
     * @return  the number of read entries to process.
     */
    public int getNumberOfReadEntriesToProcess() {
        return numberOfReadEntriesToProcess;
    }

    /**
     * Set the number of read entries to process per input file.
     * Set to <= 0 to process the entire file.
     *
     * @param numberOfReadEntriesToProcess the number of read entries to process.
     */
    public void setNumberOfReadEntriesToProcess(final int numberOfReadEntriesToProcess) {
        this.numberOfReadEntriesToProcess = numberOfReadEntriesToProcess;
    }


    /**
     * The average quality scores, one for each input file.
     * @return the list of average quality scores.
     */
    public IntList getAvgQualScores() {
        return avgQualScores;
    }

    /**
     * The names of the likely quality encoding schemes, one for each input file.
     * @return the list of likely quality encoding scheme names.
     */
    public List<String> getLikelyEncodings() {
        return likelyEncodings;
    }

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException   error parsing
     * @throws com.martiansoftware.jsap.JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        inputFilenames.clear();
        Collections.addAll(inputFilenames, jsapResult.getStringArray("input"));
        numberOfReadEntriesToProcess = jsapResult.getInt("number-of-reads");
        return this;
    }

    /**
     * Perform the conversion fasta -> compact-reads on one or more files.
     *
     * @throws java.io.IOException if the input/output files cannot be read/written
     */
    @Override
    public void execute() throws IOException {
        try {
            for (final String inputFilename : inputFilenames) {
                processingStart(inputFilename);
                final int numEntriesSampled;
                if (inputFilename.endsWith(".compact-reads")) {
                    numEntriesSampled = processCompactReadsFile(inputFilename);
                } else {
                    numEntriesSampled = processFastqFile(inputFilename);
                }
                processingEnd(numEntriesSampled);
            }
        } catch (Exception e) {
            LOG.error("Error processing", e);
        }
    }

    /**
     * Before each file, reset the state, report the filename about to be processed.
     * @param inputFilename the filename to be processed.
     */
    private void processingStart(final String inputFilename) {
        System.out.printf("Processing %s%n", inputFilename);
        numQualScoresSampled = 0;
        sumQualScores = 0;
        minQualScore = Integer.MAX_VALUE;
        maxQualScore = Integer.MIN_VALUE;
        qualityScoresFound = false;
    }

    /**
     * After each file, report the details for the file.
     * @param numEntries the number of read entries processed.
     */
    private void processingEnd(final int numEntries) {
        System.out.printf("Processed %d read entries.%n", numEntries);
        final int avgQualScore;
        if (numQualScoresSampled > 0) {
            avgQualScore = (int) (sumQualScores / numQualScoresSampled);
        } else {
            avgQualScore = 0;
        }
        final String likelyEncoding;
        if (!qualityScoresFound) {
            likelyEncoding = "fasta";
        } else if (avgQualScore <= 41) {
            likelyEncoding = "Phred";
        } else if (avgQualScore <= 83) {
            likelyEncoding = "Sanger";
        } else {
            likelyEncoding = "Illumina/Solexa";
        }

        System.out.printf("Min quality score: %d%n", minQualScore);
        System.out.printf("Max quality score: %d%n", maxQualScore);
        System.out.printf("Avg quality score: %d%n", avgQualScore);
        System.out.printf("Probable quality encoding scheme: %s%n", likelyEncoding);
        likelyEncodings.add(likelyEncoding);
        avgQualScores.add(avgQualScore);
    }

    /**
     * Process ONE compact-reads file.
     * @param inputFilename the filename to process
     * @return the number of read entries processed
     * @throws IOException error reading file
     */
    private int processCompactReadsFile(final String inputFilename) throws IOException {
        // Create directory for output file if it doesn't already exist
        int i = 0;
        for (final Reads.ReadEntry entry : new ReadsReader(inputFilename)) {
            final byte[] qualityScores = entry.getQualityScores().toByteArray();
            final boolean hasQualityScores = entry.hasQualityScores() && !ArrayUtils.isEmpty(qualityScores);
            if (hasQualityScores) {
                qualityScoresFound = true;
                for (final int qualScore : qualityScores) {
                    minQualScore = Math.min(qualScore, minQualScore);
                    maxQualScore = Math.max(qualScore, maxQualScore);
                    numQualScoresSampled++;
                    sumQualScores += Math.abs(qualScore);
                }
                if (++i == numberOfReadEntriesToProcess) {
                    break;
                }
            }
        }
        return i;
    }

    /**
     * Process ONE FASTQ file.
     * @param inputFilename the filename to process
     * @return the number of read entries processed
     * @throws IOException error reading file
     */
    private int processFastqFile(final String inputFilename) throws IOException {
        // Create directory for output file if it doesn't already exist
        int i = 0;
        for (final FastXEntry entry : new FastXReader(inputFilename)) {
            final MutableString quality = entry.getQuality();
            if (quality.length() != 0) {
                qualityScoresFound = true;
                for (final int qualScore : quality.array()) {
                    minQualScore = Math.min(qualScore, minQualScore);
                    maxQualScore = Math.max(qualScore, maxQualScore);
                    numQualScoresSampled++;
                    sumQualScores += Math.abs(qualScore);
                }
                if (++i == numberOfReadEntriesToProcess) {
                    break;
                }
            }
        }
        return i;
    }

    /**
     * Main entry point.
     * @param args command line arguments
     * @throws IOException error with IO
     * @throws JSAPException error with command line processing
     */
    public static void main(final String[] args) throws IOException, JSAPException {
        new SampleQualityScoresMode().configure(args).execute();
    }
}
