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

import it.unimi.dsi.fastutil.objects.ObjectArraySet;

/**
 *
 * User: nyasha
 * Date: Apr 8, 2010
 * Time: 12:04:01 PM
 * Calculates log2(fold change) from first group to second group (requires exactly two groups).
 */
public class Log2FoldChangeCalculator extends StatisticCalculator {

     /**
     * The natural log of the number two.
     */
    private static final double LOG_2 = Math.log(2);

        public Log2FoldChangeCalculator(final DifferentialExpressionResults results) {
            this();
            setResults(results);


        }

        public Log2FoldChangeCalculator() {
            super();
        }


        @Override
       public boolean canDo(final String[] group) {
            return group.length == 2;
        }

        @Override
       public DifferentialExpressionInfo evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                            final NormalizationMethod method, final DifferentialExpressionResults results,
                                            final DifferentialExpressionInfo info,
                                            final String... group) {
            final String groupA = group[0];
            final String groupB = group[1];
            final int foldChangeStatIndex = defineStatisticId(results, "log2-fold-change", method, group);

            final ObjectArraySet<String> samplesA = differentialExpressionCalculator.getSamples(groupA);
            final ObjectArraySet<String> samplesB = differentialExpressionCalculator.getSamples(groupB);
            double averageA = 0;
            double averageB = 0;


            for (final String sample : samplesA) {
                averageA += differentialExpressionCalculator.getNormalizedExpressionValue(sample, method, info.elementId);
            }
            for (final String sample : samplesB) {
                averageB += differentialExpressionCalculator.getNormalizedExpressionValue(sample, method, info.elementId);
            }

            final double foldChangeStatistic = log2((averageA+1) / (averageB+1));
            info.statistics.size(results.getNumberOfStatistics());
            info.statistics.set(foldChangeStatIndex, foldChangeStatistic);

            return info;
        }


    /**
     * Calculate the log2 of x
     *
     * @param x
     * @return log2(x+1)=Math.log1p(x)/Math.log(2)
     */
    private double log2(final double x) {
        return Math.log(x) / LOG_2;
    }

}
