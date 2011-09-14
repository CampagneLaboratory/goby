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
 * Converts a <a href="http://en.wikipedia.org/wiki/FASTA_format">FASTA</a>
 * or <a href="http://en.wikipedia.org/wiki/FASTQ_format">FASTQ</a> file to the compact reads
 * format. Compact reads are in the chunked protocol buffer file format described by Reads.proto.
 * Since Goby 1.7, this mode can load paired-end runs into single compact files.
 *
 * @author Fabien Campagne
 *         Date: Apr 28, 2009
 *         Time: 6:03:56 PM
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
    private static final String MODE_DESCRIPTION = "Determines the minimum and maximum quality score values " +
            "in a portion of a fastq file to help determine the likely quality encoding.";

    // The reads files to process.
    private final List<String> inputFilenames = new ArrayList<String>();

    /**
     * The CURRENT min found qual score.
     */
    private int minQualScore;

    /**
     * The CURRENT max found qual score.
     */
    private int maxQualScore;

    /**
     * The CURRENT sum of all of the quality scores, so we can perform an average.
     */
    private double sumQualScores;

    /**
     * The CURRENT number of qual scores sampled, so we can perform an average.
     */
    private int numQualScoresSampled;

    /**
     * The average qual scores for each processed file.
     */
    private final IntList avgQualScores = new IntArrayList();

    /**
     * The likely encodings for each processed file.
     */
    private final List<String> likelyEncodings = new ArrayList<String>();

    /**
     * The number of reads to to read.
     */
    private int numberOfReads = 10000;

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
     * Get the number of reads to sample.
     *
     * @return the number of reads to sample.
     */
    public int getNumberOfReads() {
        return numberOfReads;
    }

    /**
     * Set the number of reads to sample.
     *
     * @param numberOfReads the number of reads to sample.
     */
    public void setNumberOfReads(final int numberOfReads) {
        this.numberOfReads = numberOfReads;
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

    public List<String> getInputFilenames() {
        return inputFilenames;
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
        final String[] inputFilenamesArray = jsapResult.getStringArray("input");
        inputFilenames.clear();
        Collections.addAll(inputFilenames, inputFilenamesArray);
        numberOfReads = jsapResult.getInt("number-of-reads");
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

    private void processingStart(final String inputFilename) {
        System.out.printf("Processing %s%n", inputFilename);
        numQualScoresSampled = 0;
        sumQualScores = 0.0d;
        minQualScore = Integer.MAX_VALUE;
        maxQualScore = Integer.MIN_VALUE;
    }

    private void processingEnd(final int numEntries) {
        System.out.printf("Processed %d read entries.%n", numEntries);
        final int avgQualScore = (int) (sumQualScores / numQualScoresSampled);
        final String likelyEncoding;
        if (avgQualScore <= 41) {
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

    private int processCompactReadsFile(final String inputFilename) throws IOException {

        // Create directory for output file if it doesn't already exist
        int i = 0;
        for (final Reads.ReadEntry entry : new ReadsReader(inputFilename)) {
            final byte[] qualityScores = entry.getQualityScores().toByteArray();
            final boolean hasQualityScores = entry.hasQualityScores() && !ArrayUtils.isEmpty(qualityScores);
            if (hasQualityScores) {
                for (final int qualScore : qualityScores) {
                    minQualScore = Math.min(qualScore, minQualScore);
                    maxQualScore = Math.max(qualScore, maxQualScore);
                    numQualScoresSampled++;
                    sumQualScores += Math.abs(qualScore);
                }
                if (++i == numberOfReads) {
                    break;
                }
            }
        }
        return i;
    }

    private int processFastqFile(final String inputFilename) throws IOException {
        // Create directory for output file if it doesn't already exist
        int i = 0;
        for (final FastXEntry entry : new FastXReader(inputFilename)) {
            final MutableString quality = entry.getQuality();
            if (quality.length() != 0) {
                for (final int qualScore : quality.array()) {
                    minQualScore = Math.min(qualScore, minQualScore);
                    maxQualScore = Math.max(qualScore, maxQualScore);
                    numQualScoresSampled++;
                    sumQualScores += Math.abs(qualScore);
                }
                if (++i == numberOfReads) {
                    break;
                }
            }
        }
        return i;
    }

    public static void main(final String[] args) throws IOException, JSAPException {
        new SampleQualityScoresMode().configure(args).execute();
    }
}
