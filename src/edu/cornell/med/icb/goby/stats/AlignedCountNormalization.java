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

import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * Normalize by the number of alignments observed in each sample. First described in Moriatti et al 2008.
 *
 * @author Fabien Campagne
 *         Date: Mar 27, 2010
 *         Time: 4:32:40 PM
 */
public class AlignedCountNormalization extends RpkmLikeNormalizationMethod {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(AlignedCountNormalization.class);

    /**
     * {@inheritDoc}
     */
    public String getIdentifier() {
        return "aligned-count";
    }

    /**
     * {@inheritDoc}
     */
    public String getAbbreviation() {
        return "AC";
    }

    /**
     * {@inheritDoc}
     */
    public void normalize(final DifferentialExpressionCalculator calculator, final String... groups) {
        // nothing to do: deCalculator already stored the normalization factor.
        calculator.resetSumOverlapCounts();
        final ObjectSet<String> samplesToNormalize = new ObjectOpenHashSet<String>();
        for (final String group : groups) {
            samplesToNormalize.addAll(calculator.getSamples(group));
        }
        for (final String sampleId : samplesToNormalize) {
            final double normFactor = calculator.getNumAlignedInSample(sampleId);
            LOG.info(String.format("normalization denominator %g for sample %s", normFactor, sampleId));
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getDenominator(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                 final String sampleId) {                     
        return (double) differentialExpressionCalculator.getNumAlignedInSample(sampleId);
    }
}
