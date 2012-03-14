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


/**
 * Test for difference between two proportions, using the two-proportion z-test.
 * See http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cts=1331678933809&ved=0CCYQFjAA&url=http%3A%2F%2Facademic.brooklyn.cuny.edu%2Feconomic%2Ffriedman%2FTtesttwoproportionsWEB.doc&ei=0M5fT7fHGYPt0gH5iI3OBw&usg=AFQjCNHYNyIZGnOV2nN-vC0Sdb-BFuzf6Q&sig2=KXQ6_M-LCSOSq7dVw_O21A
 *
 * @author Fabien Campagne
 *         Date: 3/13/12
 *         Time: 6:51 PM
 */
public class PTestStatisticAdaptor extends AbstractMethylationAdapter {
    private static final double MAXIMUM_BOUND = 500;
    private static final long serialVersionUID = 2934190953936250446L;

    public String statName() {
        return "ptest";    // proportion t-test
    }

    LinearBinningStrategy linearBinner = new LinearBinningStrategy();

    @Override
    /**
     * Arguments must be provided in this order: Cma, Ca, Cmb, Cb.
     */
    public double calculateNoCovariate(final int... a) {
        final float cma = a[0];
        final float ca = a[1];
        final float cmb = a[2];
        final float cb = a[3];
        double pa = cma / (ca + cma);
        double pb = cmb / (cb + cmb);
        // p_hat is the pooled proportion estimate:
        double p_hat = (cma + cmb + 0.01) / (cma + cmb + ca + cb + 0.01);

        double Z = (pa - pb) /
                Math.sqrt(p_hat * (1 - p_hat) * ((1f / (cma + ca)) + (1f / (cmb + cb))));
        if (Z != Z) {
            return 0;
        } else {
            return Math.abs(Z);
        }

    }

    @Override
    /**
     *
     */
    public double calculateWithCovariate(final int covariate, final int... a) {
      final double stat= calculateNoCovariate(a);
        //  System.out.printf("stat4: stat=%g midPoint=%d %s %n", stat, midpoint, IntArrayList.wrap(a));
        return Math.min(stat, MAXIMUM_BOUND);
    }

    @Override
    /**
     * Get the upper bound on the stat adjusted for covariate.
     */
    public double getMaximumStatistic() {
        return MAXIMUM_BOUND;
    }

    @Override
    /**
     * Get the range on the stat adjusted for covariate.
     */
    public double getRange() {
        return MAXIMUM_BOUND;
    }


}
