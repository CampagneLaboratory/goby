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
 * Estimate the difference in methylation rate between sample A and B.
 * @author Fabien Campagne
 *         Date: 2/24/12
 *         Time: 2:16 PM
 */
public final class MethylationRateDifferenceStatisticAdaptor implements StatisticAdaptor {

    private static final long serialVersionUID = 3208125998662041904L;

    @Override
    public String statName() {
        return "dMR";
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


        final double totalA = cma + ca;
       final  double mra = cma /  totalA;

      final  double totalB = cmb + cb;
       final double mrb = cmb / totalB;

        return Math.abs(mra - mrb)*100.0;

    }

    @Override
    public double getMaximumStatistic() {
        return 100.0;
    }

    @Override
    public double getRange() {
        return 100.0;
    }
}
