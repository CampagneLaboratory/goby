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
import edu.cornell.med.icb.goby.algorithmic.algorithm.CoverageAnalysis;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.counts.CountsArchiveReader;
import edu.cornell.med.icb.goby.counts.CountsReader;
import edu.cornell.med.icb.goby.util.DoInParallel;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


/**
 * Write annotations corresponding to consensus peaks found in each sequence of count archives.
 */
public class CoverageMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "coverage";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Calculate various coverage statistics for a capture experiment.";
    /**
     * Basename of the counts archives to analyze.
     */
    private String[] inputBasenames;
    /**
     * Output filename where to write the stats.
     */
    private String statsOuputFilename;
    /**
     * Basename of the annotation counts archive.
     */
    private String annotationBasename;
    private int numThreads;
    private double[] percentiles = new double[]{0.9, .75, .5, .1, .01};
    private int[] depths = new int[]{5, 10, 15, 20, 30};



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
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputBasenames = jsapResult.getStringArray("input");
        statsOuputFilename = jsapResult.getString("output");
        annotationBasename = jsapResult.getString("annotation-basename");
        numThreads = jsapResult.getInt("num-threads");
        percentiles = stringToDoubles(jsapResult.getString("percentiles"));
        depths = stringToInts(jsapResult.getString("depths"));
        return this;
    }

    private double[] stringToDoubles(String percentiles) {
        String[] values = percentiles.split(",");
        DoubleArrayList result = new DoubleArrayList();
        for (String value : values) {
            result.add(Double.parseDouble(value)/100.0d);
        }
        return result.toDoubleArray();
    }

    private int[] stringToInts(String depths) {
        String[] values = depths.split(",");
        IntArrayList result = new IntArrayList();
        for (String value : values) {
            result.add(Integer.parseInt(value));
        }
        return result.toIntArray();
    }

    /**
     * Run the mode.
     */
    @Override
    public void execute() {
        final PrintWriter output;
        try {
            output = statsOuputFilename.equals("-") ? new PrintWriter(System.out) : new PrintWriter(new FileWriter(statsOuputFilename));


            DoInParallel forLoop = new DoInParallel(numThreads) {
                @Override
                public void action(DoInParallel forDataAccess, String inputBasename, int loopIndex) {
                    process(output, inputBasename);
                }
            };

            forLoop.execute(true, inputBasenames);
        } catch (IOException e) {
            System.err.println("An error occured opening the output file. ");
            System.exit(1);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void process(PrintWriter output, String basename) {
        try {
            basename = AlignmentReaderImpl.getBasename(basename);
            final AlignmentReaderImpl alignment = new AlignmentReaderImpl(basename);
            alignment.readHeader();
            alignment.close();
            final IndexedIdentifier referenceIds = alignment.getTargetIdentifiers();
            final DoubleIndexedIdentifier backwards = new DoubleIndexedIdentifier(referenceIds);

            final CountsArchiveReader archiveReader = new CountsArchiveReader(basename);
            final CountsArchiveReader annotationArchiveReader = new CountsArchiveReader(annotationBasename);

            final ObjectOpenHashSet<String> archiveIdentifiers = new ObjectOpenHashSet<String>();
            archiveIdentifiers.addAll(annotationArchiveReader.getIdentifiers());
            final CoverageAnalysis analysis = new CoverageAnalysis();

            for (int referenceIndex = 0; referenceIndex < archiveReader.getNumberOfIndices(); referenceIndex++) {
                final CountsReader reader = archiveReader.getCountReader(referenceIndex);
                // determine the corresponding chromosome in the annotation count archive:
                final String countArchiveRefId = backwards.getId(referenceIndex).toString();
                if (archiveIdentifiers.contains(countArchiveRefId)) {
                    System.out.println("Processing reference " + countArchiveRefId);
                    final CountsReader annotationReader = annotationArchiveReader.getCountReader(countArchiveRefId);

                    analysis.process(annotationReader, reader);
                } else {
                    System.out.printf("Skipping chromosome: %s%n", countArchiveRefId);
                }
            }
            final long sumDepth = analysis.getSumDepth();
            final long countDepth = analysis.getCountDepth();
            final double averageDepth = divide(sumDepth, countDepth);
            System.out.printf("Average depth= %g %n", averageDepth);
            final long sumDepthAnnot = analysis.getSumDepthAnnot();
            final long countDepthAnnot = analysis.getCountDepthAnnot();
            final double averageDepthCaptured = divide(sumDepthAnnot, countDepthAnnot);
            System.out.printf("Average depth over annotations= %g %n", averageDepthCaptured);
            analysis.estimateStatistics();

            System.out.printf("Enrichment efficiency is %2g%%%n", 100d * analysis.getEnrichmentEfficiency());
            System.out.printf("90%% of captured sites have depth>= %d%n", analysis.depthCapturedAtPercentile(.9));
            System.out.printf("75%% of captured sites have depth>= %d%n", analysis.depthCapturedAtPercentile(.75));
            System.out.printf("50%% of captured sites have depth>= %d%n", analysis.depthCapturedAtPercentile(.5));
            System.out.printf("1%% of captured sites have depth>= %d%n", analysis.depthCapturedAtPercentile(.01));
            output.printf("average-depth-captured\t%s\t%s\t%g%n", basename, "-", averageDepth);
            output.printf("average-depth\t%s\t%s\t%g%n", basename, "-", averageDepthCaptured);
            output.printf("enrichment-efficiency\t%s\t%g%%\t-%n", basename, 100d * analysis.getEnrichmentEfficiency());


            for (double percentile : percentiles) {
                output.printf("depth-captured\t%s\t%s\t%d%n", basename, Integer.toString((int) (percentile * 100)),
                        analysis.depthCapturedAtPercentile(percentile));
            }

            for (int depth : depths) {
                output.printf("percent-capture-sites-at-depth\t%s\t%s\t%d%n", basename,
                        Integer.toString((int) (100 * analysis.percentSitesCaptured(depth))),
                        depth);
            }
            output.flush();
        } catch (IOException e) {
            System.err.println("Cannot open input basename : " + basename);
            e.printStackTrace();
        }
    }

    private double sum(long[] array) {
        double sum = 0;
        int o = 0;
        for (long value : array) {

            sum += value;

        }
        return sum;
    }

    private double divide(double v1, double v2) {
        return v1 / v2;
    }

    private double sum(LongArrayList array, int offset) {
        double sum = 0;
        int o = 0;
        for (long value : array) {
            if (o >= offset) {
                sum += value;
            }
            ++o;
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
        new CoverageMode().configure(args).execute();
    }
}