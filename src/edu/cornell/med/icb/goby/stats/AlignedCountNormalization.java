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

import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.log4j.Logger;

/**
 * Normalize by the number of alignments observed in each sample.
 *
 * @author Fabien Campagne
 *         Date: Mar 27, 2010
 *         Time: 4:32:40 PM
 */
public class AlignedCountNormalization extends RpkmLikeNormalizationMethod {
    private static final Logger LOG = Logger.getLogger(AlignedCountNormalization.class);

    public String getIdentifier() {
        return "campagne-aligned-count";
    }

    public String getAbbreviation() {
        return "CAC";
    }

    public void normalize(final DifferentialExpressionCalculator calculator, final String... groups) {
        // nothing to do: deCalculator already stored the normalization factor.

        final ObjectSet<String> samplesToNormalize = new ObjectOpenHashSet<String>();
        for (final String group : groups) {
            samplesToNormalize.addAll(calculator.getSamples(group));
        }
        for (final String sampleId : samplesToNormalize) {
            final double normFactor = calculator.getNumAlignedInSample(sampleId);
            LOG.info(String.format("normalization denominator %g for sample %s", normFactor, sampleId));
        }
    }

    @Override
    public double getDenominator(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                 final String sampleId) {

        return (double) differentialExpressionCalculator.getNumAlignedInSample(sampleId);
    }
}
