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

import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: Jan 13, 2010
 *         Time: 11:35:23 AM
 */
public class AverageFisherRCalculator extends StatisticCalculator {
    public AverageFisherRCalculator(final DifferentialExpressionResults results) {
        this();
        setResults(results);
    }

    public AverageFisherRCalculator() {
        super();
    }

    /**
     * Can do as long as there is at least one group.
     */
    @Override
    public boolean canDo(final String[] group) {
        return group.length > 0;
    }

    @Override
    public DifferentialExpressionInfo evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                               final NormalizationMethod method,
                                               final DifferentialExpressionResults results,
                                               final DifferentialExpressionInfo info,
                                               final String... group) {

        final MutableString averageStatsName = new MutableString("average-Fisher-R p-value");
        final int averageStatIndex = defineStatisticId(results, averageStatsName);


        final MutableString minStatsName = new MutableString("min-Fisher-R p-value");
        final int minStatIndex = defineStatisticId(results, minStatsName);


        double averagePValue = 0;
        double minPValue = Double.MAX_VALUE;
        final IntList statIndices = results.statisticsIndexesFor("fisher-exact-R", method);

        for (final int index : statIndices) {
            final double pValue = results.getStatistic(info, index);
            averagePValue += pValue;
            minPValue = Math.min(minPValue, pValue);
        }
        averagePValue /= (double) statIndices.size();


        info.statistics.size(results.getNumberOfStatistics());
        info.statistics.set(averageStatIndex, averagePValue);
        info.statistics.set(minStatIndex, minPValue);


        return info;
    }


    public MutableString getStatisticId(final String groupId, final String modifier, final NormalizationMethod normalizationMethod) {
        return new MutableString("average " + modifier + " group " + groupId + "(" + normalizationMethod.getAbbreviation() + ")");
    }

    public MutableString getStatisticId(final String groupId, final String modifier) {
        return new MutableString("average " + modifier + " group " + groupId);
    }
}
