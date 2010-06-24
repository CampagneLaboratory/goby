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
 * Estimates sample counts, writing one new column per sample in the comparison group.
 *
 * @author Fabien Campagne
 *         Date: Jan 13, 2010
 *         Time: 11:35:23 AM
 */
public class SampleCountCalculator extends StatisticCalculator {
    public SampleCountCalculator(final DifferentialExpressionResults results) {
        this();
        setResults(results);
    }

    public SampleCountCalculator() {
        super();
    }

    /**
     * Can always do.
     */
    @Override
    public boolean canDo(final String[] group) {
        return true;
    }

    @Override
    public DifferentialExpressionInfo evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                               final NormalizationMethod method,
                                               final DifferentialExpressionResults results,
                                               final DifferentialExpressionInfo info,
                                               final String... group) {

            defineStatisticId(method, results, "RPKM", differentialExpressionCalculator.samples());
            defineStatisticId(method, results, "log2_RPKM", differentialExpressionCalculator.samples());
            defineStatisticId(null, results, "count", differentialExpressionCalculator.samples());


        for (final String sampleId : differentialExpressionCalculator.samples()) {
            final MutableString rpkmStatId = getStatisticId(sampleId, "RPKM", method);
            final int rpkmStatIndex = defineStatisticId(results, rpkmStatId);

            final MutableString log2RpkmStatId = getStatisticId(sampleId, "log2_RPKM", method);
            final int log2RpkmStatIndex = defineStatisticId(results, log2RpkmStatId);

            final MutableString countStatisticId = getStatisticId(sampleId, "count");  // counts are not normalized
            final int countStatIndex = defineStatisticId(results, countStatisticId);

                  double normalizedExpressionValue = 0;
            double rawCount = 0;

            normalizedExpressionValue += differentialExpressionCalculator.getNormalizedExpressionValue(sampleId, method, info.getElementId());
            rawCount += differentialExpressionCalculator.getOverlapCount(sampleId, info.getElementId());

            info.statistics.size(results.getNumberOfStatistics());
            info.statistics.set(rpkmStatIndex, normalizedExpressionValue);
            info.statistics.set(log2RpkmStatIndex, Math.log(normalizedExpressionValue) / Math.log(2));
            info.statistics.set(countStatIndex, rawCount);
        }

        return info;
    }

    private void defineStatisticId(NormalizationMethod method, DifferentialExpressionResults results, String statisticPrefix, String... samples) {
        for (final String sampleId : samples) {

            final MutableString statisticId = method==null? getStatisticId(sampleId, statisticPrefix) :
                    getStatisticId(sampleId, statisticPrefix, method);
            defineStatisticId(results, statisticId);
        }
    }


    public MutableString getStatisticId(final String sampleId, final String modifier, final NormalizationMethod normalizationMethod) {
        return new MutableString(modifier + " sample " + sampleId + "(" + normalizationMethod.getAbbreviation() + ")");
    }

    public MutableString getStatisticId(final String sampleId, final String modifier) {
        return new MutableString(modifier + " sample " + sampleId);
    }
}
