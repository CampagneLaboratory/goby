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
 * Estimate the Stat5 statistic, normalize by the midPoint of the sumTotal bin to remove scale effects.
 *
 * @author Fabien Campagne
 *         Date: 2/24/12
 *         Time: 2:14 PM
 */
public class Stat5StatisticAdaptor extends AbstractMethylationAdapter {
    private static final double MAXIMUM_BOUND = 20;
    private static final long serialVersionUID = 2934190953936250446L;

    public String statName() {
        return "stat5";
    }
     LinearBinningStrategy linearBinner=new LinearBinningStrategy();
    @Override
    /**
     * Arguments must be provided in this order: Cma, Ca, Cmb, Cb.
     */
    public double calculateNoCovariate(final int... a) {
        final int cma = a[0];
        final int ca = a[1];
        final int cmb = a[2];
        final int cb = a[3];
        final int diffA = cma - ca;
        final int diffB = cmb - cb;

        final int stat5;
        if (diffA * diffB > 0) {
            // change in the same direction:
            stat5 = Math.max(Math.abs(cma-cmb) , Math.abs(ca-cb));
        } else {
            //change in  opposite directions:
            stat5 = Math.abs(cma-cmb) + Math.abs(ca-cb);
        }
        return stat5;

    }

    @Override
    /**
     *
     */
    public double calculateWithCovariate(final int covariate, final int... a) {
        final double stat = calculateNoCovariate(a);
        final int binIndex = linearBinner.getBinIndex(covariate);

        // TODO would be more logical to normalize by the median of the statistic in the given bin.
        final int midpoint = linearBinner.getMidpoint(binIndex);
      //  System.out.printf("stat4: stat=%g midPoint=%d %s %n", stat, midpoint, IntArrayList.wrap(a));
        return Math.min(stat / midpoint,MAXIMUM_BOUND);
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
