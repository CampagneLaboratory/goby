/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.goby.stats.FisherExactRCalculator;
import edu.cornell.med.icb.goby.stats.StatisticsWriter;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.rosuda.JRI.Rengine;

import java.io.PrintWriter;
import java.util.Collections;
import java.util.Comparator;

/**
 * Helper class to implement the logic of discovering sequence variations in and across groups of samples.
 * Implements most of the work done by DiscoverSequenceVariantsMode
 *
 * @author Fabien Campagne
 *         Date: Sep 7, 2010
 *         Time: 2:14:38 PM
 * @see edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode
 */
public class DiscoverVariantIterateSortedAlignments
        extends IterateSortedAlignmentsListImpl {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(DiscoverVariantIterateSortedAlignments.class);

    private String[] groups;
    private ObjectArrayList<DiscoverSequenceVariantsMode.ReadIndexStats> readIndexStats;
    private int numberOfGroups;
    private int[] refCountsPerGroup;
    private int[] variantsCountPerGroup;
    private int[] distinctReadIndexCountPerGroup;
    private float[] averageVariantQualityScorePerGroup;
    private DifferentialExpressionAnalysis deAnalyzer;
    private DifferentialExpressionCalculator deCalculator;
    private int thresholdDistinctReadIndices = 10;
    private int[] readerIndexToGroupIndex;
    private int minimumVariationSupport = 3;
    private boolean fisherRInstalled;
    private int log2OddsRatioStandardErrorColumnIndex;
    private int log2OddsRatioZColumnIndex;
    private int numberOfSamples;

    private int[] refCountsPerSample;
    private int[] variantsCountPerSample;

    public void setMinimumVariationSupport(int minimumVariationSupport) {
        this.minimumVariationSupport = minimumVariationSupport;
    }

    public void setThresholdDistinctReadIndices(int thresholdDistinctReadIndices) {
        this.thresholdDistinctReadIndices = thresholdDistinctReadIndices;
    }

    public DiscoverVariantIterateSortedAlignments() {

    }

    int refIdColumnIndex;
    int positionColumnIndex;
    int log2OddsRatioColumnIndex;
    int fisherExactPValueColumnIndex;
    int observedVariationsAtPositionColumnIndex;
    StatisticsWriter statWriter;
    String[] samples;

    public void initialize(DifferentialExpressionAnalysis deAnalyzer,
                           DifferentialExpressionCalculator deCalculator,
                           String[] groups,
                           ObjectArrayList<DiscoverSequenceVariantsMode.ReadIndexStats> readIndexStats,
                           PrintWriter outWriter) {

        if (deAnalyzer.eval("within-groups") || deAnalyzer.eval("between-groups")) {
            //activate R only if we need it:
            final Rengine rEngine = GobyRengine.getInstance().getRengine();
            fisherRInstalled = rEngine != null && rEngine.isAlive();

        }
        if (deAnalyzer.eval("between-groups") && groups.length != 2) {
            System.err.println("--eval between-groups requires exactly two groups.");
            System.exit(1);
        }
        statWriter = new StatisticsWriter(outWriter);
        this.groups = groups;
        this.readIndexStats = readIndexStats;
        refIdColumnIndex = statWriter.defineColumn("referenceId");
        positionColumnIndex = statWriter.defineColumn("position");
        log2OddsRatioColumnIndex = StatisticsWriter.COLUMN_NOT_DEFINED;
        fisherExactPValueColumnIndex = -1;
        numberOfGroups = groups.length;

        samples = deCalculator.samples();
        numberOfSamples = samples.length;

        refCountsPerGroup = new int[numberOfGroups];
        variantsCountPerGroup = new int[numberOfGroups];
        distinctReadIndexCountPerGroup = new int[numberOfGroups];
        averageVariantQualityScorePerGroup = new float[numberOfGroups];

        refCountsPerSample = new int[numberOfSamples];
        variantsCountPerSample = new int[numberOfSamples];

        this.deAnalyzer = deAnalyzer;
        this.deCalculator = deCalculator;
        if (this.deAnalyzer.eval("between-groups")) {
            log2OddsRatioColumnIndex = statWriter.defineColumn("log2(Odds-ratio[%s/%s])", groups[0], groups[1]);
            log2OddsRatioStandardErrorColumnIndex = statWriter.defineColumn("log2_odds-ratio_standard_error");
            log2OddsRatioZColumnIndex = statWriter.defineColumn("log2_odds-ratio-Z");
            fisherExactPValueColumnIndex = statWriter.defineColumn("Fisher-Exact-P-value[%s/%s]", groups[0], groups[1]);
        }


        for (String group : groups) {

            statWriter.defineColumn("refProportion[%s]", group);
            statWriter.defineColumn("refCountsPerGroup[%s]", group);
            statWriter.defineColumn("varCountPerGroup[%s]", group);
            statWriter.defineColumn("distinct-read-index-count[%s]", group);
            statWriter.defineColumn("average-variant-quality-scores[%s]", group);


            if (deAnalyzer.eval("within-groups")) {
                statWriter.defineColumn("within-group-p-value[%s]", group);
            }

            statWriter.defineColumn("observed variations at position ([frequency:from/to,]+) group %s", group);

        }
       /*
        TODO uncomment to enable per sample stats.
       for (String sample : deCalculator.samples()) {

            statWriter.defineColumn("refProportion[%s]", sample);
            statWriter.defineColumn("refCountsInSample[%s]", sample);
            statWriter.defineColumn("varCountInSample[%s]", sample);
        }     */
        statWriter.writeHeader();
    }

    public void finish() {
        statWriter.close();

    }

    public void setReaderIndexToGroupIndex(int[] readerIndexToGroupIndex) {
        this.readerIndexToGroupIndex = readerIndexToGroupIndex;
    }


    public class PositionBaseInfo {
        public int readIndex;
        public int readerIndex;
        public byte qualityScore;
        public boolean matchesReference;
        public char from;
        public char to;
        public int position;
    }


    public void processPositions(int referenceIndex, int position,
                                 ObjectArrayList<IterateSortedAlignmentsListImpl.PositionBaseInfo> list) {
        int sumVariantCounts = 0;

        for (int i = 0; i < numberOfGroups; i++) {
            refCountsPerGroup[i] = 0;
            variantsCountPerGroup[i] = 0;
            distinctReadIndexCountPerGroup[i] = 0;
            averageVariantQualityScorePerGroup[i] = 0;
        }
        for (int j = 0; j < numberOfSamples; j++) {
            refCountsPerSample[j] = 0;
            variantsCountPerSample[j] = 0;

        }
        if (list != null) {
            IntSet distinctReadIndices = new IntArraySet();
            for (IterateSortedAlignmentsListImpl.PositionBaseInfo info : list) {
                //TODO make sure readerIndex matches the index of the sample name in deCalculator.samples().
                final int sampleIndex = info.readerIndex;

                final int groupIndex = readerIndexToGroupIndex[info.readerIndex];

                refCountsPerGroup[groupIndex] += info.matchesReference ? 1 : 0;
                variantsCountPerGroup[groupIndex] += info.matchesReference ? 0 : 1;

                refCountsPerSample[sampleIndex] += info.matchesReference ? 1 : 0;
                variantsCountPerSample[sampleIndex] += info.matchesReference ? 0 : 1;

                if (!info.matchesReference) {
                    sumVariantCounts += 1;
                    averageVariantQualityScorePerGroup[groupIndex] += info.qualityScore;
                    distinctReadIndices.add(info.readIndex);
                }

            }


            for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {

                averageVariantQualityScorePerGroup[groupIndex] /= variantsCountPerGroup[groupIndex];
                distinctReadIndexCountPerGroup[groupIndex] = distinctReadIndices.size();
            }

            if (distinctReadIndices.size() >= thresholdDistinctReadIndices && sumVariantCounts > minimumVariationSupport) {
                int groupIndexA = 0;
                int groupIndexB = 1;
                // Do not write statistics for positions in the start flap. The flap start is used to accumulate
                // base counts for reads that can overlap with the window under consideration.

                if (!isWithinStartFlap(referenceIndex, position)) {
                    CharSequence currentReferenceId = this.getReferenceId(referenceIndex);

                    statWriter.setValue(refIdColumnIndex, currentReferenceId);
                    statWriter.setValue(positionColumnIndex, position + 1);

                    if (deAnalyzer.eval("between-groups")) {
                        final double denominator = (double) refCountsPerGroup[groupIndexA] * (double) variantsCountPerGroup[groupIndexB];
                        double oddsRatio = denominator == 0 ? Double.NaN :
                                ((double) refCountsPerGroup[groupIndexB] * (double) variantsCountPerGroup[groupIndexA]) /
                                        denominator;
                        double logOddsRatioSE;
                        if (variantsCountPerGroup[groupIndexA] < 10 ||
                                variantsCountPerGroup[groupIndexB] < 10 ||
                                refCountsPerGroup[groupIndexA] < 10 ||
                                refCountsPerGroup[groupIndexB] < 10) {
                            // standard error estimation is unreliable when any of the counts are less than 10.
                            logOddsRatioSE = Double.NaN;
                        } else {
                            logOddsRatioSE = Math.sqrt(1d / refCountsPerGroup[groupIndexB] +
                                    1d / variantsCountPerGroup[groupIndexA] +
                                    1d / variantsCountPerGroup[groupIndexB] +
                                    1d / refCountsPerGroup[groupIndexA]);
                        }
                        double log2OddsRatio = Math.log(oddsRatio) / Math.log(2);
                        double log2OddsRatioZValue = log2OddsRatio / logOddsRatioSE;

                        double fisherP = Double.NaN;

                        boolean ok = checkCounts();
                        if (ok) {
                            fisherP = fisherRInstalled ? FisherExactRCalculator.getFisherPValue(
                                    refCountsPerGroup[groupIndexB], variantsCountPerGroup[groupIndexB],
                                    refCountsPerGroup[groupIndexA], variantsCountPerGroup[groupIndexA]) : Double.NaN;
                        } else {
                            System.err.printf("An exception was caught evaluating the Fisher Exact test P-value. Details are provided below%n" +
                                    "referenceId=%s referenceIndex=%d position=%d %n" +
                                    "refCountsPerGroup[1]=%d variantsCountPerGroup[1]=%d%n" +
                                    "refCountsPerGroup[0]=%d, variantsCountPerGroup[0]=%d",
                                    currentReferenceId, referenceIndex,
                                    position + 1,
                                    refCountsPerGroup[groupIndexB], variantsCountPerGroup[groupIndexB],
                                    refCountsPerGroup[groupIndexA], variantsCountPerGroup[groupIndexA]
                            );

                        }
                        statWriter.setValue(log2OddsRatioColumnIndex, log2OddsRatio);
                        statWriter.setValue(log2OddsRatioStandardErrorColumnIndex, logOddsRatioSE);
                        statWriter.setValue(log2OddsRatioZColumnIndex, log2OddsRatioZValue);
                        statWriter.setValue(fisherExactPValueColumnIndex, fisherP);

                    }
                    // output reference proportion, refCounts and varCounts for each sample:
                 /*
                  TODO uncomment to enable per sample stats.
                  for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
                        double refProportion = (double) refCountsPerSample[sampleIndex];
                        refProportion /= refCountsPerSample[sampleIndex] + variantsCountPerSample[sampleIndex];
                        statWriter.setValue(refProportion,
                                "refProportion[%s]", samples[sampleIndex]);
                        statWriter.setValue(refCountsPerSample[sampleIndex], "refCountsInSample[%s]", samples[sampleIndex]);
                        statWriter.setValue(variantsCountPerSample[sampleIndex], "varCountInSample[%s]", samples[sampleIndex]);

                    }  */
                    // output group information:
                    for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {
                         double refProportion = (double) refCountsPerSample[groupIndex];
                        refProportion /= refCountsPerSample[groupIndex] + variantsCountPerSample[groupIndex];
                        
                        statWriter.setValue(refProportion, "refProportion[%s]", groups[groupIndex]);
                        statWriter.setValue(refCountsPerGroup[groupIndex], "refCountsPerGroup[%s]", groups[groupIndex]);
                        statWriter.setValue(variantsCountPerGroup[groupIndex], "varCountPerGroup[%s]", groups[groupIndex]);

                        statWriter.setValue(distinctReadIndexCountPerGroup[groupIndex], "distinct-read-index-count[%s]", groups[groupIndex]);
                        statWriter.setValue(averageVariantQualityScorePerGroup[groupIndex], "average-variant-quality-scores[%s]", groups[groupIndex]);

                        if (deAnalyzer.eval("within-groups")) {
                            statWriter.setValue(estimateWithinGroupDiscoveryPalue(position, groupIndex, list, variantsCountPerGroup, refCountsPerGroup),
                                    "within-group-p-value[%s]", groups[groupIndex]);
                        }
                        summarizeVariations(statWriter, list, groupIndex);
                    }


                    statWriter.writeRecord();
                }
            }


        }
    }

    private double estimateWithinGroupDiscoveryPalue(int position, int groupIndex,
                                                     ObjectArrayList<IterateSortedAlignmentsListImpl.PositionBaseInfo> list,
                                                     int[] variantsCount, int[] refCounts) {
        double pValue = 1;
        if (!deAnalyzer.eval("within-groups")) return Double.NaN;
        int observedReferenceCount = 0;
        double expectedVariationRate = 0;
        int expectedVariationCount = 0;
        int expectedReferenceCount = 0;
        int observedVariationCount = 0;
        observedReferenceCount = refCounts[groupIndex];
        observedVariationCount = variantsCount[groupIndex];
        if (observedReferenceCount + observedVariationCount == 0) {
            // cannot call a variant if we do not observe this position in this group at all.
            return 1;
        }
        long sum = 0;
        for (IterateSortedAlignmentsListImpl.PositionBaseInfo info : list) {
            final DiscoverSequenceVariantsMode.ReadIndexStats stats = readIndexStats.get(info.readerIndex);
            final int readIndex = info.readIndex;
            if (readIndex < 1 ||
                    readIndex > stats.countVariationBases.length ||
                    readIndex > stats.countReferenceBases.length) {
                // TODO this test is a kludge: readIndex should always be within the bounds of countVariationBases or countReferenceBases
                // the fact that we need this test indicates that there is a bug in the calculation of readIndex, probably when read insertions or deletions are present.
                //              assert false : "should not fail" + String.format("readIndex =%d || readIndex >= stats.countVariationBases.length || readIndex >= stats.countReferenceBases.length",
                //                     readIndex);
                continue;
            }
            long variationBases = stats.countVariationBases[readIndex - 1];
            long referenceBases = stats.countReferenceBases[readIndex - 1];
            //  System.out.printf("readIndex=%d variationBases=%d referenceBases=%d %n",readIndexInfo.readIndex, variationBases, referenceBases);
            expectedVariationRate += variationBases;
            sum += variationBases + referenceBases;
            //   refBaseCountAnyPosition += referenceBases;
        }
        expectedVariationRate /= sum;
        //   System.out.printf("Expected variation rate: %f%n", expectedVariationRate);
        expectedVariationCount = (int) Math.round(expectedVariationRate * (double) (observedVariationCount + observedReferenceCount));
        expectedReferenceCount = (int) Math.round((1 - expectedVariationRate) * (double) (observedVariationCount + observedReferenceCount));
        if (LOG.isTraceEnabled()) {
            LOG.trace(String.format("contingency position=%d: %n" +
                    "    [  exp    obs ]%n" +
                    "ref [ %d       %d ]%n" +
                    "var [ %d       %d ] %n",
                    position,
                    expectedVariationCount, observedVariationCount,
                    expectedReferenceCount, observedReferenceCount));
        }
        pValue = fisherRInstalled ? FisherExactRCalculator.getFisherOneTailedLesserPValue(
                expectedVariationCount, observedVariationCount,
                expectedReferenceCount, observedReferenceCount
        ) : Double.NaN;
        //  System.out.printf("position=%d P-Value=%f%n", position, pValue);
        return pValue;


    }

    private boolean checkCounts() {
        boolean ok = true;
        // detect if any count is negative (that's a bug)
        for (int count : refCountsPerGroup) {

            if (count < 0) ok = false;
        }
        for (int count : variantsCountPerGroup) {
            if (count < 0) ok = false;
        }
        return ok;
    }

    private void summarizeVariations(StatisticsWriter statWriter, ObjectArrayList<IterateSortedAlignmentsListImpl.PositionBaseInfo> list,
                                     int groupIndex) {
        final Object2IntMap<MutableString> tally = new Object2IntArrayMap<MutableString>();
        tally.defaultReturnValue(0);
        for (IterateSortedAlignmentsListImpl.PositionBaseInfo info : list) {
            final int varGroupIndex = readerIndexToGroupIndex[info.readerIndex];

            if (!info.matchesReference && varGroupIndex == groupIndex) {
                MutableString variation = new MutableString();
                variation.append(info.from);
                variation.append('/');
                variation.append(info.to);
                int count = tally.getInt(variation);
                tally.put(variation, count + 1);
            }
        }
        // sort variations by decreasing tally:
        ObjectSet<MutableString> varSet = (tally.keySet());
        ObjectList<MutableString> varList = new ObjectArrayList<MutableString>();
        varList.addAll(varSet);
        Collections.sort(varList, new Comparator<MutableString>() {
            public int compare(MutableString variationA, MutableString variationB) {
                // sort by decreasing tally:
                return tally.getInt(variationB) - tally.getInt(variationA);
            }
        });

        StringBuffer buffer = new StringBuffer();
        for (MutableString variation : varList) {
            buffer.append(String.format("%d:%s,", tally.getInt(variation), variation));
        }
        statWriter.setValue(buffer.toString(),
                "observed variations at position ([frequency:from/to,]+) group %s", groups[groupIndex]);
    }
}