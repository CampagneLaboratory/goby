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
 * Estimate the delta statistic.
 *
 * @author Fabien Campagne
 *         Date: 2/24/12
 *         Time: 2:14 PM
 */
public final class Stat3StatisticAdaptor implements StatisticAdaptor {
    private static final double MAXIMUM_BOUND = 10000;
    private static final long serialVersionUID = 2934190953936250446L;

    public String statName() {
        return "stat3";
    }

    @Override
    /**
     * Arguments must be provided in this order: Cma, Ca, Cmb, Cb.
     */
    public double calculate(final int... a) {
        final int cma = a[0];
        final int ca = a[1];
        final int cmb = a[2];
        final int cb = a[3];
        final int diffA = cma - ca;
        final int diffB = cmb - cb;

        /*
       Cma=495     -90
Ca=405
Cmb=95              -90
Cb=5

and the other has:

Cma = 250
Ca = 250
Cmb= 250
Cb = 250
        */
        final int stat3;
        if (diffA * diffB > 0) {
            // change in the same direction:
            stat3 = Math.max(Math.abs(diffA) , Math.abs(diffB));
        } else {
            //change in  opposite directions:
            stat3 = Math.abs(diffA) + Math.abs(diffB);
        }

        return Math.min(stat3, MAXIMUM_BOUND);

    }

    @Override
    /**
     *
     */
    public double calculateWithCovariate(int covariate, int... a) {
        double stat = calculate(a);
        return stat / (covariate + stat + 1);
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
