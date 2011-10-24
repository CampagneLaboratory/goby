/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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
 * Estimates sample counts, writing one new column per sample. Sample ids are written without
 * any modification. This is the format expected by the Goby stats mode.
 * @author Fabien Campagne
 *         Date: Jan 13, 2010
 *         Time: 11:35:23 AM
 */
public class CountRawSampleIdsCalculator extends StatisticCalculator {
    public CountRawSampleIdsCalculator(final DifferentialExpressionResults results) {
        this();
        setResults(results);
    }

    public CountRawSampleIdsCalculator() {
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

        for (final String sampleId1 : differentialExpressionCalculator.samples()) {

            defineStatisticId(results, new MutableString(sampleId1));
        }

        for (final String sampleId : differentialExpressionCalculator.samples()) {

            final int countStatIndex = defineStatisticId(results, new MutableString(sampleId));

            double rawCount = differentialExpressionCalculator.getOverlapCount(sampleId, info.getElementId());

            info.statistics.size(results.getNumberOfStatistics());
            info.statistics.set(countStatIndex, rawCount);
        }

        return info;
    }


}
