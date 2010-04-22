/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.stats;

import edu.cornell.med.icb.goby.R.FisherExact;
import gominer.Fisher;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.inference.ChiSquareTest;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

/**
 * @author Fabien Campagne
 *         Date: Jan 12, 2010
 *         Time: 1:14:57 PM
 */
public class TestStatistics {
    @Test
    public void testFoldChange() {
        final Random randomEngine = new Random();
        final DifferentialExpressionCalculator deCalc = new DifferentialExpressionCalculator() {

            @Override
            public double getNormalizedExpressionValue(final String sample, final NormalizationMethod method, final MutableString elementId) {
                if (sample.startsWith("A")) {
                    return 2 * Math.abs(randomEngine.nextDouble());
                } else {
                    return Math.abs(randomEngine.nextDouble());
                }

                // fold change A/B = 2
            }
        };

        deCalc.defineElement("id-1");
        deCalc.defineElement("id-2");
        deCalc.defineGroup("A");
        deCalc.defineGroup("B");
        final int numReplicates = 20000;
        deCalc.reserve(2, numReplicates * 2);

        for (int i = 0; i < numReplicates; i++) {
            deCalc.associateSampleToGroup("A-" + i, "A");
            deCalc.associateSampleToGroup("B-" + i, "B");
        }
        //deCalc.associateSampleToGroup("A-", "A");
        //deCalc.associateSampleToGroup("B-1", "B");


        final DifferentialExpressionInfo info = new DifferentialExpressionInfo("id-1");
        final DifferentialExpressionResults results = new DifferentialExpressionResults();
        final FoldChangeCalculator foldChange = new FoldChangeCalculator(results);
        final NormalizationMethod normalizationMethod = new AlignedCountNormalization();
        foldChange.evaluate(deCalc, normalizationMethod, results, info, "A", "B");
        assertEquals("fold-change must be two fold", 2d, results.getStatistic(info, foldChange.statisticIds.get(0)), .1);
    }

    @Test
    public void testAverage() throws IOException {
        final Random randomEngine = new Random();
        final DifferentialExpressionCalculator deCalc = new DifferentialExpressionCalculator() {

            @Override
            public double getNormalizedExpressionValue(final String sample, final NormalizationMethod method, final MutableString elementId) {
                if (sample.startsWith("A")) {
                    return 2 * Math.abs(randomEngine.nextDouble());
                } else {
                    return Math.abs(randomEngine.nextDouble());
                }

                // fold change A/B = 2
            }
        };

        deCalc.defineElement("id-1");
        deCalc.defineElement("id-2");
        deCalc.defineGroup("A");
        deCalc.defineGroup("B");
        final int numReplicates = 20000;
        deCalc.reserve(2, numReplicates * 2);

        for (int i = 0; i < numReplicates; i++) {
            deCalc.associateSampleToGroup("A-" + i, "A");
            deCalc.associateSampleToGroup("B-" + i, "B");
        }
        //deCalc.associateSampleToGroup("A-", "A");
        //deCalc.associateSampleToGroup("B-1", "B");

        final DifferentialExpressionInfo info = new DifferentialExpressionInfo("id-1");
        final DifferentialExpressionResults results = new DifferentialExpressionResults();
        final AverageCalculator averageCalculator = new AverageCalculator(results);
        results.add(info);
        final NormalizationMethod normalizationMethod = new AlignedCountNormalization();
        averageCalculator.evaluate(deCalc, normalizationMethod, results, info, "A", "B");
        assertEquals("average A must be around 2", 1d, results.getStatistic(info, averageCalculator.getStatisticId("A", "RPKM", normalizationMethod)), .1);
        assertEquals("average B must be around 1", 0.5d, results.getStatistic(info, averageCalculator.getStatisticId("B", "RPKM", normalizationMethod)), .1);
        System.out.println(results);
        results.write(new PrintWriter("test-results/out-stats.tsv"), '\t', deCalc);
    }

    @Test
    public void testTwoStats() {
        final Random randomEngine = new Random();
        final DifferentialExpressionCalculator deCalc = new DifferentialExpressionCalculator() {

            @Override
            public double getNormalizedExpressionValue(final String sample, final NormalizationMethod method, final MutableString elementId) {
                if (sample.startsWith("A")) {
                    return 2 * Math.abs(randomEngine.nextGaussian());
                } else {
                    return Math.abs(randomEngine.nextGaussian());
                }

                // fold change A/B = 2
            }

            @Override
            public int getOverlapCount(final String sample, final MutableString elementId) {
                final NormalizationMethod normalizationMethod = new AlignedCountNormalization();
                return (int) (getNormalizedExpressionValue(sample, normalizationMethod, elementId) * 100);
            }
        };

        deCalc.defineElement("id-1");
        deCalc.defineElement("id-2");
        deCalc.defineGroup("A");
        deCalc.defineGroup("B");
        final int numReplicates = 2000;
        deCalc.reserve(2, numReplicates * 2);

        for (int i = 0; i < numReplicates; i++) {
            deCalc.associateSampleToGroup("A-" + i, "A");
            deCalc.associateSampleToGroup("B-" + i, "B");
        }

        // observe the counts to populate internal data structures:
        for (final String sampleId : deCalc.samples()) {
            final MutableString id1 = new MutableString("id-1");
            final MutableString id2 = new MutableString("id-2");
            deCalc.observe(sampleId, "id-1",
                    deCalc.getOverlapCount(sampleId, id1)
            );
            deCalc.observe(sampleId, "id-2",
                    deCalc.getOverlapCount(sampleId, id2)
            );
        }
        //deCalc.associateSampleToGroup("A-", "A");
        //deCalc.associateSampleToGroup("B-1", "B");

        final DifferentialExpressionInfo info = new DifferentialExpressionInfo("id-1");
        final DifferentialExpressionResults results = new DifferentialExpressionResults();
        final FoldChangeCalculator foldChange = new FoldChangeCalculator(results);
        final TTestCalculator tTest = new TTestCalculator(results);
        final FisherExactTestCalculator fisher = new FisherExactTestCalculator(results);
        final NormalizationMethod normalizationMethod = new AlignedCountNormalization();
        foldChange.evaluate(deCalc, normalizationMethod, results, info, "A", "B");
        tTest.evaluate(deCalc, normalizationMethod, results, info, "A", "B");
        fisher.evaluate(deCalc, normalizationMethod, results, info, "A", "B");
        assertEquals("fold-change must be two fold", 2d, results.getStatistic(info, foldChange.statisticIds.get(0)), .1);
        assertTrue("T-test must be significant", results.getStatistic(info, tTest.statisticIds.get(0)) < 0.01);
        assertTrue("fisher test must not be significant", results.getStatistic(info, fisher.statisticIds.get(0)) > 0.05);
    }

    @Test
    public void testFisher() throws MathException {
        final DifferentialExpressionCalculator deCalc = new DifferentialExpressionCalculator();
        final int numReplicates = 2;
        deCalc.defineElement("id-1");
        deCalc.defineElement("id-2");
        deCalc.defineGroup("A");
        deCalc.defineGroup("B");
        deCalc.reserve(2, numReplicates * 2);

        for (int i = 1; i <= numReplicates; i++) {
            deCalc.associateSampleToGroup("A-" + i, "A");
            deCalc.associateSampleToGroup("B-" + i, "B");
        }

        /**
         * Encode the following table in two genes:
         Fisher's Exact Test
         http://www.langsrud.com/fisher.htm
         ------------------------------------------
         TABLE = [ 10 , 20 , 30 , 40 ]
         Left   : p-value = 0.2533310713617698
         Right  : p-value = 0.8676419647894328
         2-Tail : p-value = 0.5044757698516504
         ------------------------------------------
         */
        deCalc.observe("A-1", "id-1", 7);
        deCalc.observe("A-2", "id-1", 3);         // 7+3 = 10
        deCalc.observe("B-1", "id-1", 15);
        deCalc.observe("B-2", "id-1", 5);         // 15+5 =20

        deCalc.observe("A-1", "id-2", 15);
        deCalc.observe("A-2", "id-2", 15);        // 15+15=30
        deCalc.observe("B-1", "id-2", 20);
        deCalc.observe("B-2", "id-2", 20);        // 20+20=40

        final DifferentialExpressionInfo info = new DifferentialExpressionInfo("id-1");
        final DifferentialExpressionResults results = new DifferentialExpressionResults();
        final FisherExactTestCalculator fisher = new FisherExactTestCalculator(results);
        final NormalizationMethod normalizationMethod = new AlignedCountNormalization();
        fisher.evaluate(deCalc, normalizationMethod, results, info, "A", "B");
        assertEquals("fisher test equal expected result", 0.5044757698516504, results.getStatistic(info, fisher.statisticIds.get(0)), 0.001);


        final Fisher fisherTest = new Fisher();
        final int totalCountInA = 1700;
        final int totalCountInB = 170; // equal total in each group
        final int sumCountInA = 90;
        final int sumCountInB = 45; // half the counts in sample B

        fisherTest.fisher(totalCountInA, sumCountInA, totalCountInA + totalCountInB, sumCountInA + sumCountInB);

        final double pValue = fisherTest.getTwotail();
        final double proportionTotalA = divide(totalCountInA, (totalCountInA + totalCountInB));
        final double proportionTotalB = divide(totalCountInB, (totalCountInA + totalCountInB));
        final ChiSquareTest chisquare = new ChiSquareTestImpl();
        final double nGroups = 2;
        final double[] expected = {divide(sumCountInA + sumCountInB, nGroups) * proportionTotalA * nGroups,
                divide(sumCountInA + sumCountInB, nGroups) * proportionTotalB * nGroups};
        final long[] observed = {sumCountInA, sumCountInB};
        final double chiPValue = Math.abs(chisquare.chiSquareTest(expected, observed));

        assertTrue("pValue: " + chiPValue, chiPValue < 0.001);
// The Fisher implementation we are using return 1 for the above. This is wrong. Compare to the chi-square result
// (results should be comparable since the counts in each cell are large)
//         assertTrue("pValue: " + pValue, pValue < 0.001);
    }

    @Test
    public void testFisherExact() throws MathException {
        final DifferentialExpressionCalculator deCalc = new DifferentialExpressionCalculator();
        final int numReplicates = 2;
        deCalc.defineElement("id-1");
        deCalc.defineElement("id-2");
        deCalc.defineGroup("A");
        deCalc.defineGroup("B");
        deCalc.reserve(2, numReplicates * 2);

        for (int i = 1; i <= numReplicates; i++) {
            deCalc.associateSampleToGroup("A-" + i, "A");
            deCalc.associateSampleToGroup("B-" + i, "B");
        }

        /**
         * Encode the following table in two genes:
         Fisher's Exact Test
         http://www.langsrud.com/fisher.htm
         ------------------------------------------
         TABLE = [ 10 , 20 , 30 , 40 ]
         Left   : p-value = 0.2533310713617698
         Right  : p-value = 0.8676419647894328
         2-Tail : p-value = 0.5044757698516504
         ------------------------------------------
         */
        deCalc.observe("A-1", "id-1", 7);
        deCalc.observe("A-2", "id-1", 3);         // 7+3 = 10
        deCalc.observe("B-1", "id-1", 15);
        deCalc.observe("B-2", "id-1", 5);         // 15+5 =20

        deCalc.observe("A-1", "id-2", 15);
        deCalc.observe("A-2", "id-2", 15);        // 15+15=30
        deCalc.observe("B-1", "id-2", 20);
        deCalc.observe("B-2", "id-2", 20);        // 20+20=40

        final DifferentialExpressionInfo info = new DifferentialExpressionInfo("id-1");
        final DifferentialExpressionResults results = new DifferentialExpressionResults();

        final FisherExactRCalculator fisher = new FisherExactRCalculator(results);
        if (fisher.installed()) {
            final NormalizationMethod normalizationMethod = new AlignedCountNormalization();
            fisher.evaluate(deCalc, normalizationMethod, results, info, "A", "B");
            assertEquals("fisher test equal expected result", 0.5044757698516504, results.getStatistic(info, fisher.statisticIds.get(0)), 0.001);


            final int totalCountInA = 1700;
            final int totalCountInB = 170; // equal total in each group
            final int sumCountInA = 90;
            final int sumCountInB = 45; // half the counts in sample B

            final int sumCountNotInA = totalCountInA - sumCountInA;
            final int sumCountNotInB = totalCountInB - sumCountInB;

            final FisherExact.Result result =
                    FisherExact.fexact(sumCountInA, sumCountNotInA, sumCountInB, sumCountNotInB);
            final double pValue = result.getPValue();

            final double proportionTotalA = divide(totalCountInA, (totalCountInA + totalCountInB));
            final double proportionTotalB = divide(totalCountInB, (totalCountInA + totalCountInB));
            final ChiSquareTest chisquare = new ChiSquareTestImpl();
            final double nGroups = 2;
            final double[] expected = {divide(sumCountInA + sumCountInB, nGroups) * proportionTotalA * nGroups,
                    divide(sumCountInA + sumCountInB, nGroups) * proportionTotalB * nGroups};
            final long[] observed = {sumCountInA, sumCountInB};
            final double chiPValue = Math.abs(chisquare.chiSquareTest(expected, observed));

            assertTrue("pValue: " + chiPValue, chiPValue < 0.001);
// The Fisher implementation we are using return 1 for the above. This is wrong. Compare to
// the chi-square result
// (results should be comparable since the counts in each cell are large)
            assertTrue("pValue: " + pValue, pValue < 0.001);
        }
    }

    private double divide(final int a, final int b) {
        return ((double) a) / ((double) b);
    }

    private double divide(final int a, final double b) {
        return ((double) a) / b;
    }

    @Test
    public void testChiSquare() throws MathException {

        final DifferentialExpressionCalculator deCalc = new DifferentialExpressionCalculator();
        final int numReplicates = 2;
        deCalc.defineElement("id-1");
        deCalc.defineElement("id-2");
        deCalc.defineGroup("A");
        deCalc.defineGroup("B");
        deCalc.reserve(2, numReplicates * 2);

        for (int i = 1; i <= numReplicates; i++) {
            deCalc.associateSampleToGroup("A-" + i, "A");
            deCalc.associateSampleToGroup("B-" + i, "B");
        }

        deCalc.observe("A-1", "id-1", 7);
        deCalc.observe("A-2", "id-1", 3);         // 7+3 = 10
        deCalc.observe("B-1", "id-1", 15);
        deCalc.observe("B-2", "id-1", 5);         // 15+5 =20

        deCalc.observe("A-1", "id-2", 15);
        deCalc.observe("A-2", "id-2", 15);        // 15+15=30
        deCalc.observe("B-1", "id-2", 20);
        deCalc.observe("B-2", "id-2", 20);        // 20+20=40


        final DifferentialExpressionInfo info = new DifferentialExpressionInfo("id-1");
        final DifferentialExpressionResults results = new DifferentialExpressionResults();

        final ChiSquareTestCalculator calc = new ChiSquareTestCalculator(results);
        final NormalizationMethod normalizationMethod = new AlignedCountNormalization();
        calc.evaluate(deCalc, normalizationMethod, results, info, "A", "B");
        assertEquals("chi square test equal expected result", 0.456056540250256, results.getStatistic(info, calc.statisticIds.get(0)), 0.001);

        final ChiSquareTest chisquare = new ChiSquareTestImpl();
        final double[] expected = {30, 12};
        final long[] observed = {0, 100};
        final double chiPValue = chisquare.chiSquareTest(expected, observed);

        assertTrue("pValue: " + chiPValue, chiPValue < 0.001);
// The Fisher implementation we are using return 1 for the above. This is wrong. Compare to the chi-square result
// (results should be comparable since the counts in each cell are large)
//         assertTrue("pValue: " + pValue, pValue < 0.001);
    }

    @Test

    public void testChiSquareZeroCount() throws MathException {

        final DifferentialExpressionCalculator deCalc = new DifferentialExpressionCalculator();
        final int numReplicates = 2;
        deCalc.defineElement("id-1");
        deCalc.defineElement("id-2");
        deCalc.defineGroup("A");
        deCalc.defineGroup("B");
        deCalc.reserve(2, numReplicates * 2);

        for (int i = 1; i <= numReplicates; i++) {
            deCalc.associateSampleToGroup("A-" + i, "A");
            deCalc.associateSampleToGroup("B-" + i, "B");
        }

        deCalc.observe("A-1", "id-1", 0);                   // ZERO counts should yield NaN p-value
        deCalc.observe("A-2", "id-1", 0);         // 0+0 = 0
        deCalc.observe("B-1", "id-1", 15);
        deCalc.observe("B-2", "id-1", 5);         // 15+5 =20

        deCalc.observe("A-1", "id-2", 15);
        deCalc.observe("A-2", "id-2", 15);        // 15+15=30
        deCalc.observe("B-1", "id-2", 20);
        deCalc.observe("B-2", "id-2", 20);        // 20+20=40


        final DifferentialExpressionInfo info = new DifferentialExpressionInfo("id-1");
        final DifferentialExpressionResults results = new DifferentialExpressionResults();

        final ChiSquareTestCalculator calc = new ChiSquareTestCalculator(results);
        final NormalizationMethod normalizationMethod = new AlignedCountNormalization();
        calc.evaluate(deCalc, normalizationMethod, results, info, "A", "B");
        assertTrue("chi square test result must be NaN (zero count)", Double.isNaN(results.getStatistic(info, calc.statisticIds.get(0))));

    }


    @Test
    public void testFDR() {
        final Random randomEngine = new Random();
        randomEngine.setSeed(1013);
        final BonferroniAdjustment bonferroni = new BonferroniAdjustment();
        final BenjaminiHochbergAdjustment fdr = new BenjaminiHochbergAdjustment();
        final DifferentialExpressionResults list = new DifferentialExpressionResults();
        final String statId = "t-test-P-value";
        list.declareStatistic(statId);
        final int statIndex = list.getStatisticIndex(statId);
        final int numObservations = 100000;
        final double proportionOfNaN = .1;
        for (int i = 0; i < numObservations; i++) {
            final DifferentialExpressionInfo info = new DifferentialExpressionInfo("element-" + i);
            info.statistics.size(list.getNumberOfStatistics());
            final double random1 = randomEngine.nextDouble();
            final double random2 = randomEngine.nextDouble();

            info.statistics.set(statIndex, random1 < proportionOfNaN ? Double.NaN : random2);
            list.add(info);
        }
        final String secondPValueId = "another-p-value";
        list.declareStatistic(secondPValueId);
        final int statIndex2 = list.getStatisticIndex(secondPValueId);
        for (final DifferentialExpressionInfo info : list) {
            info.statistics.size(list.getNumberOfStatistics());
            info.statistics.set(statIndex2, randomEngine.nextDouble());
        }
        final NormalizationMethod normalizationMethod = new AlignedCountNormalization();
        bonferroni.adjust(list, normalizationMethod, statId, secondPValueId);
        fdr.adjust(list, normalizationMethod, statId, secondPValueId);
        final int index1 = list.getStatisticIndex("t-test-P-value-BH-FDR-q-value");
        final int index2 = list.getStatisticIndex(secondPValueId + "-BH-FDR-q-value");

        final double significanceThreshold = 0.05;
        int numRejectedHypothesesTest1 = 0;
        int numRejectedHypothesesTest2 = 0;
        for (final DifferentialExpressionInfo info : list) {

            final boolean test1 = info.statistics.getDouble(index1) > significanceThreshold;
            if (!test1) {
                //   System.out.println("info:" + info);
                numRejectedHypothesesTest1++;
            }
            final boolean test2 = info.statistics.getDouble(index2) > significanceThreshold;
            if (!test2) {
                // System.out.println("info:" + info);
                numRejectedHypothesesTest2++;
            }


        }
        assertTrue("No q-value should be significant after FDR adjustment", numRejectedHypothesesTest1 < significanceThreshold * numObservations);
        assertTrue("No q-value should be significant after FDR adjustment", numRejectedHypothesesTest2 < significanceThreshold * numObservations);

        //      System.out.println("list.adjusted: " + list);


        final double[] p = {
                2.354054e-07, 2.101590e-05, 2.576842e-05, 9.814783e-05, 1.052610e-04
                , 1.241481e-04, 1.325988e-04, 1.568503e-04, 2.254557e-04, 3.795380e-04
                , 6.114943e-04, 1.613954e-03, 3.302430e-03, 3.538342e-03, 5.236997e-03
                , 6.831909e-03, 7.059226e-03, 8.805129e-03, 9.401040e-03, 1.129798e-02
                , 2.115017e-02, 4.922736e-02, 6.053298e-02, 6.262239e-02, 7.395153e-02
                , 8.281103e-02, 8.633331e-02, 1.190654e-01, 1.890796e-01, 2.058494e-01
                , 2.209214e-01, 2.856000e-01, 3.048895e-01, 4.660682e-01, 4.830809e-01
                , 4.921755e-01, 5.319453e-01, 5.751550e-01, 5.783195e-01, 6.185894e-01
                , 6.363620e-01, 6.448587e-01, 6.558414e-01, 6.885884e-01, 7.189864e-01
                , 8.179539e-01, 8.274487e-01, 8.971300e-01, 9.118680e-01, 9.437890e-01};

        final double[] adjusted_R = {
                1.177027e-05, 4.294736e-04, 4.294736e-04, 9.471343e-04, 9.471343e-04
                , 9.471343e-04, 9.471343e-04, 9.803146e-04, 1.252532e-03, 1.897690e-03
                , 2.779520e-03, 6.724807e-03, 1.263693e-02, 1.263693e-02, 1.745666e-02
                , 2.076243e-02, 2.076243e-02, 2.445869e-02, 2.473958e-02, 2.824495e-02
                , 5.035754e-02, 1.118804e-01, 1.304633e-01, 1.304633e-01, 1.479031e-01
                , 1.592520e-01, 1.598765e-01, 2.126168e-01, 3.259994e-01, 3.430823e-01
                , 3.563248e-01, 4.462501e-01, 4.619538e-01, 6.835770e-01, 6.835770e-01
                , 6.835770e-01, 7.188450e-01, 7.414352e-01, 7.414352e-01, 7.626063e-01
                , 7.626063e-01, 7.626063e-01, 7.626063e-01, 7.824868e-01, 7.988737e-01
                , 8.802645e-01, 8.802645e-01, 9.304775e-01, 9.304775e-01, 9.437890e-01};

        final double[] adjusted_R_nocummin = {
                1.177027e-05, 5.253976e-04, 4.294736e-04, 1.226848e-03, 1.052610e-03
                , 1.034567e-03, 9.471343e-04, 9.803146e-04, 1.252532e-03, 1.897690e-03
                , 2.779520e-03, 6.724807e-03, 1.270165e-02, 1.263693e-02, 1.745666e-02
                , 2.134972e-02, 2.076243e-02, 2.445869e-02, 2.473958e-02, 2.824495e-02
                , 5.035754e-02, 1.118804e-01, 1.315934e-01, 1.304633e-01, 1.479031e-01
                , 1.592520e-01, 1.598765e-01, 2.126168e-01, 3.259994e-01, 3.430823e-01
                , 3.563248e-01, 4.462501e-01, 4.619538e-01, 6.853944e-01, 6.901156e-01
                , 6.835770e-01, 7.188450e-01, 7.567830e-01, 7.414352e-01, 7.732368e-01
                , 7.760512e-01, 7.676890e-01, 7.626063e-01, 7.824868e-01, 7.988737e-01
                , 8.890803e-01, 8.802645e-01, 9.345104e-01, 9.304775e-01, 9.437890e-01
        };


        final int n = p.length;
        for (int rank = p.length; rank >= 1; rank--) {
            final int index = rank - 1;
            assertEquals("rank: " + rank, adjusted_R_nocummin[index], p[index] * (((double) n) / (double) rank), 0.01);
        }

        final DifferentialExpressionResults list2 = new DifferentialExpressionResults();
        list2.declareStatistic("p-value");
        int i = 0;
        for (final double pValue : p) {
            final DifferentialExpressionInfo info = new DifferentialExpressionInfo(String.valueOf(i++));
            info.statistics.add(pValue);
            list2.add(info);
        }

        final DifferentialExpressionResults list3 = fdr.adjust(list2, "p-value");
        System.out.println("list3:" + list3);
        final int index = list3.getStatisticIndex("p-value-BH-FDR-q-value");
        for (final DifferentialExpressionInfo infoAdjusted : list3) {
            final int elementIndex = Integer.parseInt(infoAdjusted.getElementId().toString());
            assertEquals("adjusted p-values must match for i=" + infoAdjusted.getElementId(),
                    adjusted_R[elementIndex], infoAdjusted.statistics.get(index), 0.01);
        }
    }
}
