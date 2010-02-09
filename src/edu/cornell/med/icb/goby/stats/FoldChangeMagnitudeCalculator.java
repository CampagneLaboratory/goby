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
 * Calculates fold change magnitude from first group to second group (requires exactly two groups).
 * The fold change magnitude is defined as max(foldChange, 1d / foldChange).
 *
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 7:06:31 PM
 */
public class FoldChangeMagnitudeCalculator extends StatisticCalculator {

    public FoldChangeMagnitudeCalculator(final DifferentialExpressionResults results) {
        this();
        setResults(results);


    }

    public FoldChangeMagnitudeCalculator() {
        super("fold-change-magnitude");
    }


    @Override
    boolean canDo(final String[] group) {
        return group.length == 2;
    }

    @Override
    DifferentialExpressionInfo evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                        final DifferentialExpressionResults results,
                                        final DifferentialExpressionInfo info,
                                        final String... group) {
        final String groupA = group[0];
        final String groupB = group[1];

        final ObjectArraySet<String> samplesA = differentialExpressionCalculator.getSamples(groupA);
        final ObjectArraySet<String> samplesB = differentialExpressionCalculator.getSamples(groupB);
        double averageA = 0;
        double averageB = 0;


        for (final String sample : samplesA) {
            averageA += differentialExpressionCalculator.getRPKM(sample, info.elementId);
        }
        for (final String sample : samplesB) {
            averageB += differentialExpressionCalculator.getRPKM(sample, info.elementId);
        }

        final double foldChange = averageA / averageB;
        final double foldChangeMagnitudeStatistic = Math.max(foldChange, 1d / foldChange);

        info.statistics.size(results.getNumberOfStatistics());
        info.statistics.set(results.getStatisticIndex(statisticId), foldChangeMagnitudeStatistic);

        return info;
    }


}
