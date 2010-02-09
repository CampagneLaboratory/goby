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
        super(null);
    }

    /**
     * Can do as long as there is at least one group.
      */
    @Override
    boolean canDo(final String[] group) {
        return group.length > 0;
    }

    @Override
    DifferentialExpressionInfo evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                        final DifferentialExpressionResults results,
                                        final DifferentialExpressionInfo info,
                                        final String... group) {
        for (final String groupId : group) {
            final MutableString rpkmStatId = getStatisticId(groupId, "RPKM");
            if (!results.isStatisticDefined(rpkmStatId)) {
                results.declareStatistic(rpkmStatId);
            }
            final MutableString countStatisticId = getStatisticId(groupId, "count");
            if (!results.isStatisticDefined(countStatisticId)) {
                results.declareStatistic(countStatisticId);
            }

            // calculate the average over the group:
            final ObjectArraySet<String> samplesA = differentialExpressionCalculator.getSamples(groupId);

            double averageRPKM = 0;
            double averageCount = 0;

            for (final String sample : samplesA) {
                averageRPKM += differentialExpressionCalculator.getRPKM(sample, info.elementId);
                averageCount += differentialExpressionCalculator.getOverlapCount(sample, info.elementId);
            }
            averageRPKM /= (double) samplesA.size();
            averageCount /= (double) samplesA.size();

            info.statistics.size(results.getNumberOfStatistics());
            info.statistics.set(results.getStatisticIndex(rpkmStatId), averageRPKM);
            info.statistics.set(results.getStatisticIndex(countStatisticId), averageCount);
        }

        return info;
    }

    public MutableString getStatisticId(final String groupId, final String modifier) {
        return new MutableString("average " + modifier + " group " + groupId);
    }
}
