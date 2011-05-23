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
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.FileWriter;


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

        return this;
    }

    /**
     * Run the mode.
     */
    @Override
    public void execute() {
        PrintWriter output = null;
        try {
            output = statsOuputFilename.equals("-") ? new PrintWriter(System.out) : new PrintWriter(new FileWriter(statsOuputFilename));

        } catch (IOException e) {
            System.err.println("An error occured opening the output file. ");
            System.exit(1);
        }
        for (String basename : inputBasenames) {
            try {
                final AlignmentReaderImpl alignment = new AlignmentReaderImpl(basename);
                alignment.readHeader();
                alignment.close();
                final IndexedIdentifier referenceIds = alignment.getTargetIdentifiers();
                final DoubleIndexedIdentifier backwards = new DoubleIndexedIdentifier(referenceIds);

                CountsArchiveReader archiveReader = new CountsArchiveReader(basename);
                CountsArchiveReader annotationArchiveReader = new CountsArchiveReader(annotationBasename);

                ObjectOpenHashSet<String> archiveIdentifiers = new ObjectOpenHashSet<String>();
                archiveIdentifiers.addAll(annotationArchiveReader.getIdentifiers());
                CoverageAnalysis analysis = new CoverageAnalysis();

                for (int referenceIndex = 0; referenceIndex < archiveReader.getNumberOfIndices(); referenceIndex++) {
                    CountsReader reader = archiveReader.getCountReader(referenceIndex);
                    // determine the corresponding chromosome in the annotation count archive:
                    String countArchiveRefId = backwards.getId(referenceIndex).toString();
                    if (archiveIdentifiers.contains(countArchiveRefId)) {
                        System.out.println("Processing reference " + countArchiveRefId);
                        CountsReader annotationReader = annotationArchiveReader.getCountReader(countArchiveRefId);

                        analysis.process(annotationReader, reader);
                    } else {
                        System.out.printf("Skipping chromosome: %s%n", countArchiveRefId);
                    }
                }
                long sumDepth = analysis.getSumDepth();
                long countDepth = analysis.getCountDepth();
                final double averageDepth = divide(sumDepth, countDepth);
                System.out.printf("Average depth= %g %n", averageDepth);
                long sumDepthAnnot = analysis.getSumDepthAnnot();
                long countDepthAnnot = analysis.getCountDepthAnnot();
                final double averageDepthCaptured = divide(sumDepthAnnot, countDepthAnnot);
                System.out.printf("Average depth over annotations= %g %n", averageDepthCaptured);
                analysis.estimateStatistics();

                System.out.printf("Enrichment efficiency cumulative is %2g%%%n", 100d * analysis.getEnrichmentEfficiency());
                System.out.printf("90%% of captured sites have depth>= %d%n", analysis.depthCapturedAtPercentile(.9));
                System.out.printf("75%% of captured sites have depth>= %d%n", analysis.depthCapturedAtPercentile(.75));
                System.out.printf("50%% of captured sites have depth>= %d%n", analysis.depthCapturedAtPercentile(.5));
                System.out.printf("1%% of captured sites have depth>= %d%n", analysis.depthCapturedAtPercentile(.01));
                output.printf("average-depth-captured\t%s\t%s\t%g%n", basename, "-", averageDepth);
                output.printf("average-depth\t%s\t%s\t%g%n", basename, "-", averageDepthCaptured);

                for (double percentile : new double[]{0.9, .75, .5, 1}) {
                    output.printf("depth-captured\t%s\t%s\t%d%n", basename, Integer.toString((int) (percentile * 100)),
                            analysis.depthCapturedAtPercentile(percentile));
                }
                output.flush();
            } catch (IOException e) {
                System.err.println("Cannot open input basename : " + basename);
                e.printStackTrace();
            }
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