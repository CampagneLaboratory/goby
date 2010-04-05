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
 * Implement an RPKM calculation with a custom denominator (effective number of reads per lane).
 * @author Fabien Campagne
 *         Date: Mar 28, 2010
 *         Time: 12:48:14 PM
 */
public abstract class RpkmLikeNormalizationMethod implements NormalizationMethod {

    abstract double getDenominator(DifferentialExpressionCalculator differentialExpressionCalculator, String sampleId);

    /**
     * {@inheritDoc}
     */
    public double getNormalizedExpressionValue(final DifferentialExpressionCalculator deCalc,
                                               final String sampleId,
                                               final MutableString elementId) {
        final int count = deCalc.getOverlapCount(sampleId, elementId);
        final int elementLength = deCalc.getElementLength(elementId);

        final double normalizationFactor = getDenominator(deCalc, sampleId); // in reads
        return (double) count / ((double) elementLength / 1000.0d) / (normalizationFactor / 1E6d);
    }
}
