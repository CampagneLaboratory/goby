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

/**
 * Implement the Bonferroni multiple testing correction.
 *
 * @author Fabien Campagne
 *         Date: Jan 13, 2010
 *         Time: 11:18:00 AM
 * @see <a href="http://en.wikipedia.org/wiki/Bonferroni_correction">Wikipedia</a>
 */
public class BonferroniAdjustment extends FDRAdjustment {
    @Override
    public DifferentialExpressionResults adjust(final DifferentialExpressionResults list, final String statisticId) {
        final MutableString statistic = new MutableString(statisticId);

        final String adjustedStatisticId = statisticId + "-Bonferroni-adjusted";
        list.declareStatistic(adjustedStatisticId);
        final int adjustedStatisticIndex = list.getStatisticIndex(new MutableString(adjustedStatisticId));

        final int statisticIndex = list.getStatisticIndex(statistic);
        final double listSize = getListSize(list, statisticIndex);


        for (final DifferentialExpressionInfo info : list) {

            final double pValue = info.statistics.get(statisticIndex);
            double adjustedPValue = pValue * listSize;

            if (adjustedPValue > 1) {
                adjustedPValue = 1;
            }
            info.statistics.size(list.getNumberOfStatistics());
            info.statistics.set(adjustedStatisticIndex, adjustedPValue > 1 ? 1 : adjustedPValue);

        }
        return list;

    }


}
