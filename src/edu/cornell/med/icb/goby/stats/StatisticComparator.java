/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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

import it.unimi.dsi.lang.MutableString;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Compare two {@link edu.cornell.med.icb.goby.stats.DifferentialExpressionInfo} elements
 * by the specified statistic.
 *
 * @author Fabien Campagne
 *         Date: Jan 12, 2010
 *         Time: 6:43:53 PM
 */
public class StatisticComparator implements Comparator<DifferentialExpressionInfo>, Serializable {
    /**
     * Used for serialization.
     */
    private static final long serialVersionUID = 1402700497385045251L;
    private final int statisticIndex;

    public StatisticComparator(final DifferentialExpressionResults list,
                               final MutableString statisticId) {

        this.statisticIndex = list.getStatisticIndex(statisticId);
        assert statisticIndex != -1 : String.format("could not find statistic %s to adjust.", statisticId);
    }

    public int compare(final DifferentialExpressionInfo info1,
                       final DifferentialExpressionInfo info2) {
        final double statisticValue1 = info1.statistics.getDouble(statisticIndex);
        final double statisticValue2 = info2.statistics.getDouble(statisticIndex);
        return Double.compare(statisticValue1, statisticValue2);
    }
}
