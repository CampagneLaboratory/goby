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

import it.unimi.dsi.lang.MutableString;

import java.util.Comparator;

/**
 * Compare two DE info elements by the specified statistic.
 * @author Fabien Campagne
 *         Date: Jan 12, 2010
 *         Time: 6:43:53 PM
 */
public class StatisticComparator implements Comparator<DifferentialExpressionInfo> {
   private final  int statisticIndex;

    public StatisticComparator(DifferentialExpressionResults list, MutableString statisticId) {
        this.statisticIndex = list.getStatisticIndex(statisticId);
    }

    public int compare(DifferentialExpressionInfo info1, DifferentialExpressionInfo info2) {
        final double statisticValue1 = info1.statistics.get(statisticIndex);
        final double statisticValue2 = info2.statistics.get(statisticIndex);
        if (statisticValue1<statisticValue2) return -1;
        if (statisticValue1>statisticValue2) return 1;
        else return 0;

        //final double test = statisticValue1 - statisticValue2;
        //return (int) test;

    }
}
