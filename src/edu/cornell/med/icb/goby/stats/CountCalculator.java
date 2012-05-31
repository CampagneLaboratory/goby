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

import java.util.Map;

import org.apache.commons.lang.StringUtils;

/**
 * Estimates sample counts, writing one new column per sample in the comparison group.
 * This is similar to SampleCountCalculator but ONLY outputs counts, no RPKM and log2 RPKM.
 * Useful for creating DESeq-compliant count tables. Also, this version includes
 * the GROUP NAME in brackets in the column headers as the last part of the column header.
 *
 * @author Fabien Campagne
 *         Date: Jan 13, 2010
 *         Time: 11:35:23 AM
 */
public class CountCalculator extends StatisticCalculator {
    public CountCalculator(final DifferentialExpressionResults results) {
        this();
        setResults(results);
    }

    public CountCalculator() {
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

        defineStatisticId(null, results, "count", differentialExpressionCalculator.getSampleToGroupMap(),
                differentialExpressionCalculator.samples());

        for (final String sampleId : differentialExpressionCalculator.samples()) {

            final MutableString countStatisticId =  getStatisticId(sampleId, "count", null,
                            differentialExpressionCalculator.getSampleToGroupMap().get(sampleId)); // counts are not normalized
            final int countStatIndex = defineStatisticId(results, countStatisticId);

            double normalizedExpressionValue = 0;
            double rawCount = 0;

            normalizedExpressionValue += differentialExpressionCalculator.getNormalizedExpressionValue(sampleId, method, info.getElementId());
            rawCount += differentialExpressionCalculator.getOverlapCount(sampleId, info.getElementId());

            info.statistics.size(results.getNumberOfStatistics());
            info.statistics.set(countStatIndex, rawCount);
        }

        return info;
    }

    //
    private void defineStatisticId(
            final NormalizationMethod method, final DifferentialExpressionResults results,
            final String statisticPrefix, Map<String, String> sampleToGroupNameMap, final String... samples) {
        for (final String sampleId : samples) {
            final MutableString statisticId =
                    getStatisticId(sampleId, statisticPrefix, method, sampleToGroupNameMap.get(sampleId));
            defineStatisticId(results, statisticId);
        }
    }


    public MutableString getStatisticId(
            final String sampleId, final String modifier, final NormalizationMethod normalizationMethod, final String groupName) {
        MutableString sb = new MutableString();
        sb.append(modifier).append(" sample ").append(sampleId);
        if (normalizationMethod != null) {
            sb.append(" (").append(normalizationMethod.getAbbreviation()).append(')');
        }
        if (groupName != null) {
            sb.append(" [").append(groupName).append(']');
        }
        return sb;
    }
}
