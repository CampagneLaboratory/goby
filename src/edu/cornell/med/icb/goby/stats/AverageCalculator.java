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

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: Jan 13, 2010
 *         Time: 11:35:23 AM
 */
public class AverageCalculator extends StatisticCalculator {
    public AverageCalculator(final DifferentialExpressionResults results) {
        this();
        setResults(results);
    }

    public AverageCalculator() {
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
        // group statistic by type rather than by DE 'group':
        
        for (final String groupId : group) {
            final MutableString rpkmStatId = getStatisticId(groupId, "RPKM", method);
            defineStatisticId(results, rpkmStatId);


        }
        for (final String groupId : group) {


            final MutableString log2RpkmStatId = getStatisticId(groupId, "log2_RPKM", method);
            defineStatisticId(results, log2RpkmStatId);


        }
        for (final String groupId : group) {


            final MutableString countStatisticId = getStatisticId(groupId, "count");
            defineStatisticId(results, countStatisticId);
        }


        for (final String groupId : group) {
            final MutableString rpkmStatId = getStatisticId(groupId, "RPKM", method);
            final int rpkmStatIndex = defineStatisticId(results, rpkmStatId);

            final MutableString log2RpkmStatId = getStatisticId(groupId, "log2_RPKM", method);
            final int log2RpkmStatIndex = defineStatisticId(results, log2RpkmStatId);

            final MutableString countStatisticId = getStatisticId(groupId, "count");
            final int countStatIndex = defineStatisticId(results, countStatisticId);

            // calculate the average over the group:
            final ObjectArraySet<String> samplesA = differentialExpressionCalculator.getSamples(groupId);

            double averageNormalizedExpressionValue = 0;
            double averageCount = 0;

            for (final String sample : samplesA) {
                averageNormalizedExpressionValue += differentialExpressionCalculator.getNormalizedExpressionValue(sample, method, info.getElementId());
                averageCount += differentialExpressionCalculator.getOverlapCount(sample, info.getElementId());
            }
            averageNormalizedExpressionValue /= (double) samplesA.size();
            averageCount /= (double) samplesA.size();

            info.statistics.size(results.getNumberOfStatistics());
            info.statistics.set(rpkmStatIndex, averageNormalizedExpressionValue);
            info.statistics.set(log2RpkmStatIndex, Math.log(averageNormalizedExpressionValue) / Math.log(2));
            info.statistics.set(countStatIndex, averageCount);
        }

        return info;
    }


    public MutableString getStatisticId(final String groupId, final String modifier, final NormalizationMethod normalizationMethod) {
        return new MutableString("average " + modifier + " group " + groupId + "(" + normalizationMethod.getAbbreviation() + ")");
    }

    public MutableString getStatisticId(final String groupId, final String modifier) {
        return new MutableString("average " + modifier + " group " + groupId);
    }
}
