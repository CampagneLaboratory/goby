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
 * Do not use, in the source tree for junit compatibility.
 *
 * @author Fabien Campagne
 *         Date: 2/24/12
 *         Time: 2:14 PM
 */
public final class BuggyDeltaStatisticAdaptor implements StatisticAdaptor {
    private static final double MAXIMUM_BOUND = 10000.0;
    private static final long serialVersionUID = 2934190953936250446L;

    @Override
    public String statName() {
        return "delta(deprecated,bug)";
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

        final int maxA;
        final int maxB;
        final int minA;
        final int minB;

        if (cma > ca) {
            maxA = cma;
            minA = ca;
            maxB = cmb;
            minB = cb;
        } else {
            maxA = ca;
            maxB = cb;
            minA = cma;
            minB = cmb;
        }
        final int diffA = maxA - minA;
        final int diffB = maxB - minB;
        return Math.min(Math.abs(diffA - diffB), MAXIMUM_BOUND);

    }

    @Override
    public double getMaximumStatistic() {
        return MAXIMUM_BOUND;
    }

    @Override
    public double getRange() {
        return MAXIMUM_BOUND;
    }
}
