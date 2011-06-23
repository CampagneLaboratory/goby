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

import java.util.Collections;

/**
 * Compute Benjamini Hochberg adjusted q-values.
 *
 * @author Fabien Campagne
 *         Date: Jan 12, 2010
 *         Time: 6:39:50 PM
 * @see <a href="http://en.wikipedia.org/wiki/False_discovery_rate">Wikipedia</a>
 */
public class BenjaminiHochbergAdjustment extends FDRAdjustment {

    @Override
    public DifferentialExpressionResults adjust(final DifferentialExpressionResults list, final String statisticId) {
        final MutableString statistic = new MutableString(statisticId);
        final String adjustedStatisticId = statisticId + "-BH-FDR-q-value";
        list.declareStatistic(adjustedStatisticId);
        final int adjustedStatisticIndex = list.getStatisticIndex(new MutableString(adjustedStatisticId));
        // sort differentially expressed elements by increasing statistic (typically a P-value):

        Collections.sort(list, new StatisticComparator(list, statistic));

        final int statisticIndex = list.getStatisticIndex(statistic);
        final double listSize = getListSize(list);

        ///int rank = 1;
        double cummin = 1;

        int size = list.size();
        int completeSize = size+ignoredElementsAboveThreshold;
        int rank = completeSize;
        for (int index = size; index >= 1; --index) {
            //   for (DifferentialExpressionInfo info : list) {
            final DifferentialExpressionInfo info = list.get(index - 1);
            final double pValue = info.statistics.get(statisticIndex);
            double adjustedPValue = 1;
            if (pValue == pValue) {

                // pValue is a number.
                final double adjustment = completeSize / (double) rank;
                adjustedPValue = pValue * adjustment;
                if (adjustedPValue < cummin) {
                    // keep track of the smallest adjusted P-value seen so far (traversing from large to small P-values):
                    cummin = adjustedPValue;

                } else {
                    // if the current adjustedPvalue would be larger than a previously seen P-value, keep the previous one
                    // so that the null hypothesis is also rejected for the current element, at any significance threshold.
                    // This behaviour mimics the logic implemented in the R p.adjust( "BH") method with the cummin function.
                    adjustedPValue = cummin;
                }
            } else {
                // we just encountered a NaN p-value, reset cummin..
                cummin = 1;
            }

            if (adjustedPValue > 1) {
                adjustedPValue = 1;
            }
            /*   System.out.println(String.format("Adjusting p-value %g by listSize=%g rank=%d factor=%g => %g", pValue,
                     listSize, rank, adjustment, adjustedPValue));
            */
            info.statistics.size(list.getNumberOfStatistics());
            info.statistics.set(adjustedStatisticIndex, adjustedPValue);
            --rank;
        }
        return list;
    }

    /**
     * Set the number of elements that were not stored in list, because their P-value was already above threshold.
     *
     * @param ignoredElementsAboveThreshold the number of elements not in list.
     */
    public void setNumberOfElementsAboveThreshold(int ignoredElementsAboveThreshold) {
        //To change body of created methods use File | Settings | File Templates.
    }
}
