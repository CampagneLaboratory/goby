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

import java.util.Collections;

/**
 * Implement the Bonferroni multiple testing correction.
 * @see <a href="http://en.wikipedia.org/wiki/Bonferroni_correction">Wikipedia</a>
 * @author Fabien Campagne
 *         Date: Jan 13, 2010
 *         Time: 11:18:00 AM
 */
public class BonferroniAdjustment extends FDRAdjustment {
    public DifferentialExpressionResults adjust(DifferentialExpressionResults list, String statisticId) {
        MutableString statistic = new MutableString(statisticId);
        final String adjustedStatisticId = statisticId + "-Bonferroni-adjusted";
        list.declareStatistic(adjustedStatisticId);
        int adjustedStatisticIndex = list.getStatisticIndex(new MutableString(adjustedStatisticId));
        // sort differentially expressed elements by increasing statistic (typically a P-value):
        Collections.sort(list, new StatisticComparator(list, statistic));

        double listSize = list.size();
        final int statisticIndex = list.getStatisticIndex(statistic);
        for (DifferentialExpressionInfo info : list) {

            final double pValue = info.statistics.get(statisticIndex);
            final double adjustedPValue = pValue * listSize;
            info.statistics.size(list.getNumberOfStatistics());
            info.statistics.set(adjustedStatisticIndex, adjustedPValue);

        }
        return list;

    }
}
