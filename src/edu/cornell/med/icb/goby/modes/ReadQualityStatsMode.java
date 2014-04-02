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

import com.google.protobuf.ByteString;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.math.random.MersenneTwister;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 * Evaluate statistics for read qualities.
 *
 * @author Fabien Campagne
 */
public class ReadQualityStatsMode extends AbstractGobyMode {

    private static final org.apache.log4j.Logger LOG = org.apache.log4j.Logger.getLogger(ReadQualityStatsMode.class);
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "read-quality-stats";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Calculate statistics for quality scores in a compact reads file.";

    /**
     * The output file.
     */
    private File outputFile;

    /**
     * The basename of the compact read files.
     */
    private final List<File> inputFiles = new LinkedList<File>();
    private double sampleFraction = 0.01;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    enum OutputFormat {
        TSV,
    }

    private OutputFormat outputFormat = OutputFormat.TSV;

    /**
     * Number of samples ignored because of sampleFraction value < 1.0d.
     */
    private long numberOfSkippedReads;

    /**
     * Number of samples observed because of sampleFraction value < 1.0d. If sampleFraction == 1.0d
     * this should be the total number of samples in the file.
     */
    private long numberOfObservedReads;

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
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        final File[] inputFilesArray = jsapResult.getFileArray("input");
        inputFiles.addAll(Arrays.asList(inputFilesArray));
        outputFile = jsapResult.getFile("output");
        outputFormat = OutputFormat.valueOf(jsapResult.getString("format").toUpperCase());
        sampleFraction = jsapResult.getDouble("sample-fraction");
        return this;
    }

    /**
     * The percentage of reads to process. 0.01 means 1% of reads,
     * 1.0 means 100% of reads. The default of 0.01 should work fine
     * for most files but if you are dealing with a very small file
     * you should set this to 1.0.
     *
     * @return The percentage of reads to process
     */
    public double getSampleFraction() {
        return sampleFraction;
    }

    /**
     * The precentage of reads to process. 0.01 means 1% of reads,
     * 1.0 means 100% of reads. The default of 0.01 should work fine
     * for most files but if you are dealing with a very small file
     * you should set this to 1.0.
     *
     * @param sampleFraction The percentage of reads to process
     */
    public void setSampleFraction(final double sampleFraction) {
        this.sampleFraction = sampleFraction;
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
     * Get the output file.
     *
     * @return the output file
     */
    public File getOutputFile() {
        return outputFile;
    }

    /**
     * Set the output file.
     *
     * @param outputFile the output file
     */
    public void setOutputFile(final File outputFile) {
        this.outputFile = outputFile;
    }

    /**
     * Number of reads skipped because of sampleFraction value < 1.0d.
     * @return the number of skipped samples.
     */
    public long getNumberOfSkippedReads() {
        return numberOfSkippedReads;
    }

    /**
     * Number of reads observed because of sampleFraction value < 1.0d. If sampleFraction == 1.0d
     * this should be the total number of reads in the file.
     * @return the number of observed samples
     */
    public long getNumberOfObservedReads() {
        return numberOfObservedReads;
    }

    /**
     * Display the alignments as text files.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        PrintStream writer = null;
        try {
            numberOfSkippedReads = 0;
            numberOfObservedReads = 0;
            writer = outputFile == null ? System.out
                    : new PrintStream(new FileOutputStream(outputFile));
            final Int2ObjectMap<ReadQualityStats> qualityStats = new Int2ObjectOpenHashMap<ReadQualityStats>();
            final ProgressLogger progress = new ProgressLogger(LOG);
            progress.start();
            final MersenneTwister random = new MersenneTwister(37);
            final boolean doSample = sampleFraction < 1.0d;
            writer.println("basename\treadIndex\t25%-percentile\tmedian\taverageQuality\t75%-percentile");
            for (final File filename : inputFiles) {
                final ReadsReader reader = new ReadsReader(filename);
                // we do getName to remove any path, that is not taken out by ReadsReader.getBasename().
                final String basename = FilenameUtils.getName(ReadsReader.getBasename(filename.toString()));
                for (final Reads.ReadEntry entry : reader) {
                    if (!doSample || random.nextDouble() < sampleFraction) {
                        final ByteString qualityScores = entry.getQualityScores();
                        final int size = qualityScores.size();
                        for (int readIndex = 0; readIndex < size; readIndex++) {
                            final byte code = qualityScores.byteAt(readIndex);
                            ReadQualityStats stats = qualityStats.get(readIndex);
                            if (stats == null) {
                                stats = new ReadQualityStats(1.0d);
                                qualityStats.put(readIndex, stats);
                                stats.readIndex = readIndex;
                            }
                            stats.observe(code);
                        }
                        numberOfObservedReads++;
                    } else {
                        numberOfSkippedReads++;
                    }
                    progress.lightUpdate();
                }

                for (final ReadQualityStats stat : qualityStats.values()) {
                    if (!stat.sampleIsEmpty()) {
                        stat.evaluatePercentiles();
                        writer.printf("%s\t%d\t%d\t%d\t%f\t%d%n",
                                basename,
                                stat.readIndex,
                                stat.percentile(25),
                                stat.percentile(50),
                                stat.averageQuality / stat.observedCount,
                                stat.percentile(75));
                    }
                }
            }
            progress.updateAndDisplay();
            progress.stop();

        } finally {
            if (writer != System.out) {
                IOUtils.closeQuietly(writer);
            }
        }
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
        new ReadQualityStatsMode().configure(args).execute();
    }

    private static class ReadQualityStats {
        private final MersenneTwister random = new MersenneTwister();
        private int readIndex;
        // byte minValue;
        // byte maxValue;
        private final ByteArrayList sample = new ByteArrayList();
        private double averageQuality;
        private int observedCount;
        private final double sampleFraction;
        private final boolean doSample;

        private ReadQualityStats(final double sampleFraction) {
            this.sampleFraction = sampleFraction;
            doSample = sampleFraction < 1.0;
        }

        void observe(final byte b) {
            averageQuality += b;
            observedCount += 1;
            if (!doSample || random.nextDouble() < sampleFraction) {
                sample.add(b);
            }
        }

        public void evaluatePercentiles() {
            Collections.sort(sample);
        }

        public byte percentile(final double percent) {
            final int index = (int) ((double) sample.size() * percent / 100d);
            return sample.get(index);
        }

        public boolean sampleIsEmpty() {
            return sample.isEmpty();
        }
    }
}
