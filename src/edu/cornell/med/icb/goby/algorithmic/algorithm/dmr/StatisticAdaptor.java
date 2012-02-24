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

import java.io.Serializable;

/**
 * @author Fabien Campagne
 *         Date: 2/24/12
 *         Time: 2:09 PM
 */
public interface StatisticAdaptor extends Serializable {
    /**
     *
     * Return the name of the statistic estimated by this adaptor.
     * @return
     */
    public String statName();
    /**
     * Estimate a statistic given some integers.
     *
     * @param a values used to calculate the statistic. The meaning of the values and their order is implementation dependent.
     * @return the calculated statistic.
     */

    public double calculate(int... a);

    /**
     * Return the maximum value that the statistic can attain.
     * @return a maximum bound on statistic.
     */
    double getMaximumStatistic();

    /**
     * Get the range of this statistic. The difference between the maximum statistic that can be calculated and the
     * minimum value.
     * @return range of the statistic.
     */
    double getRange();
}
