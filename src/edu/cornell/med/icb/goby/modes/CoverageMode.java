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
import edu.cornell.med.icb.goby.counts.*;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;

import java.io.*;

import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;


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
        for (String basename : inputBasenames) {
            try {
                final AlignmentReaderImpl alignment = new AlignmentReaderImpl(basename);
                alignment.readHeader();
                alignment.close();
                final IndexedIdentifier referenceIds = alignment.getTargetIdentifiers();
                final DoubleIndexedIdentifier backwards = new DoubleIndexedIdentifier(referenceIds);

                CountsArchiveReader archiveReader = new CountsArchiveReader(basename);
                CountsArchiveReader annotationArchiveReader = new CountsArchiveReader(annotationBasename);
                /* An array where each element is the number of bases observed for which exactly i reads match span the base. The index of the array is i.

                 */
                LongArrayList depthTallyInAnnotation = new LongArrayList(10000);
                LongArrayList depthTallyOutsideAnnotation = new LongArrayList(10000);

                ObjectOpenHashSet<String> archiveIdentifiers = new ObjectOpenHashSet<String>();
                archiveIdentifiers.addAll(annotationArchiveReader.getIdentifiers());
                // sum of depth when depth is not zero:
                long sumDepth = 0;
                int countDepth = 0;
                long countAllBases = 0;
                for (int referenceIndex = 0; referenceIndex < archiveReader.getNumberOfIndices(); referenceIndex++) {
                    CountsReader reader = archiveReader.getCountReader(referenceIndex);
                    // determine the corresponding chromosome in the annotation count archive:
                    String countArchiveRefId = backwards.getId(referenceIndex).toString();
                    if (archiveIdentifiers.contains(countArchiveRefId)) {
                        System.out.println(countArchiveRefId);
                        CountsReader annotationReader = annotationArchiveReader.getCountReader(countArchiveRefId);
                        AnyTransitionCountsIterator orIterator = new AnyTransitionCountsIterator(reader, annotationReader);

                        while (orIterator.hasNextTransition()) {
                            orIterator.nextTransition();
                            int annotationCount = orIterator.getCount(1);
                            int readerCount = orIterator.getCount(0);
                            int position = orIterator.getPosition();
                            int length = orIterator.getLength();
                            int end = position + length;

                            boolean inAnnotation = annotationCount == 1;
                            LongArrayList update = inAnnotation ? depthTallyInAnnotation : depthTallyOutsideAnnotation;
                            int depth = readerCount;
                            if (depth != 0) {
                                sumDepth += depth;
                                countDepth++;
                            }
                            grow(depthTallyInAnnotation, depth);
                            grow(depthTallyOutsideAnnotation, depth);

                            update.set(depth, update.get(depth) + length);
                            countAllBases += length;

                        }
                    } else {
                        System.out.printf("Skipping chromosome: %s%n", countArchiveRefId);
                    }
                }
                System.out.printf("Average depth= %g %n", divide(sumDepth, countDepth));
                System.out.printf("Enrichment efficiency is %2g%%%n", 100d * divide(sum(depthTallyInAnnotation, 1), sum(depthTallyOutsideAnnotation, 1) + sum(depthTallyInAnnotation, 1)));
                //      System.out.println("capturedDepths: " + depthTallyInAnnotation);
                //    System.out.println("notcapture Depths: " + depthTallyOutsideAnnotation);
                /*double[] fractionOfBasesCovered = new double[depthTallyInAnnotation.size()];
                for (int i = 0; i < fractionOfBasesCovered.length; i++) {
                    fractionOfBasesCovered[i] = ((double) depthTallyInAnnotation.get(i)) /
                            ((double) (depthTallyOutsideAnnotation.get(i) + depthTallyInAnnotation.get(i)));
                }
                */
                //     System.out.println("fraction of bases covered at depths: " + DoubleArrayList.wrap(fractionOfBasesCovered));


                final int length = depthTallyInAnnotation.size();
                long[] cumulativeCaptured = new long[length];
                long[] cumulativeNotCaptured = new long[length];

                for (int depth = length - 2; depth >= 0; --depth) {
                    cumulativeCaptured[depth] = cumulativeCaptured[depth + 1] + depthTallyInAnnotation.get(depth);
                    cumulativeNotCaptured[depth] = cumulativeNotCaptured[depth + 1] + depthTallyOutsideAnnotation.get(depth);
                }

                System.out.printf("Enrichment efficiency cumulative is %2g%%%n", 100d * divide(cumulativeCaptured[0], countAllBases));
                double[] fractionOfBasesCoveredCumulative = new double[length];
                for (int i = 0; i < length; i++) {
                    fractionOfBasesCoveredCumulative[i] = ((double) cumulativeCaptured[i]) /
                            ((double) (cumulativeCaptured[i] + cumulativeNotCaptured[i]));
                }
                //      System.out.println("fraction of bases covered at >=depths : " + DoubleArrayList.wrap(fractionOfBasesCoveredCumulative));

                for (int depth = 0; depth < 100; depth++) {
                    //System.out.printf("%d %g %n", depth, divide(cumulativeCaptured[depth],cumulativeCaptured[0]) * 100);
                    System.out.printf("%d %g %g %n", depth,
                            divide(cumulativeCaptured[depth], countAllBases) * 100,
                            divide(cumulativeNotCaptured[depth], countAllBases) * 100);
                }
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

    private void grow(LongArrayList update, int depth) {
        for (int k = update.size(); k <= depth; k++) {
            update.add(0);
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
        new CoverageMode().configure(args).execute();
    }
}