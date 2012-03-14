/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author Fabien Campagne
 *         Date: 2/20/12
 *         Time: 4:26 PM
 */
public class TestDensityEstimator {
    class value {
        int context;
        int cma;
        int ca;
        int cmb;
        int cb;

        value(int context, int cma, int ca, int cmb, int cb) {
            this.context = context;
            this.cma = cma;
            this.ca = ca;
            this.cmb = cmb;
            this.cb = cb;
        }
    }

    @Test
    public void testObserve() throws Exception {
        final DensityEstimator estimator = new DensityEstimator(2, new BuggyDeltaStatisticAdaptor());
        estimator.setBinningStrategy(new LinearBinningStrategy());
        final value[] observations = {
                new value(0, 5, 13, 4, 12),
                new value(0, 0, 1, 1, 2),
                new value(0, 0, 1, 0, 3),
        };
        for (final value obs : observations) {
            estimator.observe(obs.context, obs.cma, obs.ca, obs.cmb, obs.cb);
        }
        assertEquals("", 2l, estimator.getCumulativeCount(0, 0));
        assertEquals("", 3l, estimator.getCumulativeCount(2, 0));
    }

    @Test
    public void testObserveLargeSumTotals() throws Exception {
        final DensityEstimator estimator = new DensityEstimator(2, new BuggyDeltaStatisticAdaptor());
        final value[] observations = {
                new value(0, 500, 130, 4, 12), // sumTotal= 646
                new value(0, 0, 1000, 1, 2),   // sumTotal= 1003
                new value(0, 0, 10000, 0, 3),  // sumTotal= 10003
        };
        for (final value obs : observations) {
            estimator.observe(obs.context, obs.cma, obs.ca, obs.cmb, obs.cb);
        }
        assertEquals("", 1l, estimator.getCumulativeCount( 378,646));
        assertEquals("", 1l, estimator.getCumulativeCount( 378,647));  // sumTotal is binned, and 647 goes in same bin as 646.
        assertEquals("", 0l, estimator.getCumulativeCount(  377,647));  // cumulative count before observation is zero
        assertEquals("", 1l, estimator.getCumulativeCount( 999,1003));
        assertEquals("", 1l, estimator.getCumulativeCount(9997,10003));
    }

    @Test
    public void testFastIndex() {
        LinearBinningStrategy binning = new LinearBinningStrategy();
        int[] sumTotalValues = {
                0, 1, 3, 50, 99, 101, 502, 1050, 2000, 1000010
        };
        for (int sumTotal = 0; sumTotal < 10000; sumTotal++) {
            int theIndex = binning.getBinIndex(sumTotal);
            int theIndexFast = binning.getTheIndex(sumTotal);
            assertEquals(String.format("theIndex=%d must match theIndexFast=%d for sumTotal=%d.", theIndex, theIndexFast, sumTotal), theIndex, theIndexFast);
        }

    }


    @Test
    public void testScaling() {
        StatisticAdaptor adapter = new AbstractMethylationAdapter() {
            private static final long serialVersionUID = 8506302569020149425L;

            @Override
            public String statName() {
                return "dummy";  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            public double calculateNoCovariate(final int... a) {
                return a[0];
            }

            @Override
            public double calculateWithCovariate(int covariate, int... a) {
                return calculateNoCovariate(a);
            }

            @Override
            public double getMaximumStatistic() {
                return 10;
            }

            @Override
            public double getRange() {
                return 10;
            }
        };
        DensityEstimator estimator = new DensityEstimator(1, adapter);
        int[] statistics = {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3};

        for (int stat : statistics) {
            estimator.observe(0, stat);
        }
        assertEquals(0, estimator.getCumulativeCount(0.0, 1));
        assertEquals(6, estimator.getCumulativeCount(1.0, 1));
        assertEquals(10, estimator.getCumulativeCount(2.0, 2));
        assertEquals(13, estimator.getCumulativeCount(3.0, 3));
        assertEquals(13, estimator.getCumulativeCount(4.0, 4));
        assertEquals(13, estimator.getCumulativeCount(5.0, 5));
    }

    @Test
    public void testDeltas() {
        DeltaStatisticAdaptor adaptor = new DeltaStatisticAdaptor();
        assertEquals(180.0, adaptor.calculateNoCovariate(495, 405, 95, 5), 0.1);
        assertEquals(0.0, adaptor.calculateNoCovariate(250, 250, 250, 250), 0.1);
    }

    @Test
    public void testdMR() {
        StatisticAdaptor adaptor = new MethylationRateDifferenceStatisticAdaptor();
        assertEquals(40.0, adaptor.calculateNoCovariate(495, 405, 95, 5), .1);
        assertEquals(0.0, adaptor.calculateNoCovariate(250, 250, 250, 250), 0.1);
    }

    @Test
    public void testStat3() {
        StatisticAdaptor adaptor = new Stat3StatisticAdaptor();
        assertEquals(90.0, adaptor.calculateNoCovariate(495, 405, 95, 5), .1);
        assertEquals(180.0, adaptor.calculateNoCovariate(405, 495, 95, 5), .1);
        assertEquals(0.0, adaptor.calculateNoCovariate(250, 250, 250, 250), 0.1);
        assertEquals(2.0, adaptor.calculateNoCovariate(251, 250, 252, 250), 0.1);
    }

    @Test
    public void testStat4() {
        StatisticAdaptor adaptor = new Stat4StatisticAdaptor();
        assertEquals(1.8, adaptor.calculateWithCovariate(30, 495, 405, 95, 5), .1);
        assertEquals(0.12, adaptor.calculateWithCovariate(1100, 405, 495, 95, 5), .1);
        assertEquals(0, adaptor.calculateWithCovariate(1100, 250, 250, 250, 250), 0.001);
        assertEquals(0.0013333333333333333, adaptor.calculateWithCovariate(1100, 251, 250, 252, 250), 0.001);
    }

    @Test
    public void testStat5() {
        StatisticAdaptor adaptor = new Stat5StatisticAdaptor();
        assertEquals(400.0 / 1500.0, adaptor.calculateWithCovariate(1000, 495, 405, 95, 5), .1);   // 55->90
        assertEquals(800.0 / 1500.0, adaptor.calculateWithCovariate(1100, 405, 495, 95, 5), .1);  //45->95
        assertEquals(0.0, adaptor.calculateWithCovariate(1100, 250, 250, 250, 250), 0.001);
        assertEquals(2.0 / 1500.0, adaptor.calculateWithCovariate(1100, 251, 250, 252, 250), 0.001);
      // note that stat5 is sensitive to differences in  coverage in each group, ie., it will detect when Cma+Ca differ from Cmb+Cb.
        // We introduced ptest to remove this problem.
        assertEquals(0.9514285714285714, adaptor.calculateWithCovariate(10010, 10000, 0, 10, 0), 0.001);
    }


    @Test
    public void testPTest() {
        StatisticAdaptor adaptor = new PTestStatisticAdaptor();
        assertEquals(7.715497971503937, adaptor.calculateWithCovariate(1000, 495, 405, 95, 5), .1);   // 55->90
        assertEquals(9.486833263709139, adaptor.calculateWithCovariate(1100, 405, 495, 95, 5), .1);  //45->95
        assertEquals(0, adaptor.calculateWithCovariate(1100, 250, 250, 250, 250), 0.001);
        assertEquals(0.031481125356480476, adaptor.calculateWithCovariate(1100, 251, 250, 252, 250), 0.001);
        // ptest has no problem with differences in coverage between groups:
        assertEquals(0, adaptor.calculateWithCovariate(10010, 10000, 0, 10, 0), 0.001);
    }
    @Test
       public void testFisherExact() {
           StatisticAdaptor adaptor = new FisherExactTestAdaptor();
           assertEquals(16.789604574725345, adaptor.calculateNoCovariate( 495, 405, 95, 5), .1);   // 55->90
           assertEquals(16.789604574725345, adaptor.calculateWithCovariate(1000, 495, 405, 95, 5), .1);   // 55->90
           assertEquals(23.839098600218055, adaptor.calculateWithCovariate(1100, 405, 495, 95, 5), .1);  //45->95
           assertEquals(0.0d, adaptor.calculateWithCovariate(1100, 250, 250, 250, 250), 0.001d);
           assertEquals(0.00d, adaptor.calculateWithCovariate(1100, 251, 250, 252, 250), 0.001d);
       }
    @Test
    public void testSmallLog10Binning() {
        SmallAndLog10BinningStrategy binner = new SmallAndLog10BinningStrategy();
        double covariates[] = {0, 1, 10, 40, 99, 100, 150, 1200, 10001, 100001, 1000001, 10000001, 100000001, 1000000001};
        int expectedBinIndex[] = {0, 0, 0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8};
        assertEquals(covariates.length, expectedBinIndex.length);
        for (int i = 0; i < covariates.length; i++) {
            final int binIndex = binner.getBinIndex(covariates[i]);
            assertEquals(expectedBinIndex[i], binIndex);
            assertTrue(covariates[i] >= binner.getLowerBound(binIndex));
            assertTrue(covariates[i] <= binner.getUpperBound(binIndex));
        }

    }
    @Test
    public void testFastLog10Binning() {
        FastSmallAndLog10BinningStrategy binner = new FastSmallAndLog10BinningStrategy();
        double covariates[] = {0, 1, 10, 40, 99, 100, 150, 1200, 10001, 100001, 1000001, 10000001, 100000001, 1000000001};
        int expectedBinIndex[] = {0, 0, 0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8};
        assertEquals(covariates.length, expectedBinIndex.length);
        for (int i = 0; i < covariates.length; i++) {
            final int binIndex = binner.getBinIndex(covariates[i]);
            assertEquals(expectedBinIndex[i], binIndex);
            assertTrue(covariates[i] >= binner.getLowerBound(binIndex));
            assertTrue(covariates[i] <= binner.getUpperBound(binIndex));
        }

    }
}
