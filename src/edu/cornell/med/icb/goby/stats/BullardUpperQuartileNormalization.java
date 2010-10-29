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

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.util.Collections;

/**
 * Implements the upper quartile normalization method described by Bullard et al in BMC
 * Bioinformatics 2010, 11:94 doi:10.1186/1471-2105-11-94.  This method multiples counts
 * in each sample i by a factor d_i. d_i is taken to be the 1/q_i with q_i the 75% quantile,
 * or count value larger than 75% of counts in sample i. Elements that have zero counts in
 * all samples are not considered to estimate the 75% quartile.
 *
 * @author Fabien Campagne
 *         Date: Mar 27, 2010
 *         Time: 3:03:54 PM
 */
public class BullardUpperQuartileNormalization extends RpkmLikeNormalizationMethod {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(BullardUpperQuartileNormalization.class);

    private final double percentile = .75d;
    private final Object2DoubleMap<String> normalizationFactors =
            new Object2DoubleOpenHashMap<String>();

    /**
     * {@inheritDoc}
     */
    public String getIdentifier() {
        return "bullard-upper-quartile";
    }

    /**
     * {@inheritDoc}
     */
    public String getAbbreviation() {
        return "BUQ";
    }

    /**
     * {@inheritDoc}
     */
    public void normalize(final DifferentialExpressionCalculator calculator, final String... groups) {
        calculator.resetSumOverlapCounts();
        final ObjectSet<String> samplesToNormalize = new ObjectOpenHashSet<String>();
        for (final String group : groups) {
            samplesToNormalize.addAll(calculator.getSamples(group));
        }
        // determine set of elements with reads in at least one sample (lane): elementIdsToKeep
        final ObjectSet<MutableString> allElementIds = calculator.getElementIds();
        final ObjectSet<MutableString> elementIdsToKeep = new ObjectOpenHashSet<MutableString>();
        for (final String sample : samplesToNormalize) {
            for (final MutableString elementId : allElementIds) {
                final int overlapCount = calculator.getOverlapCount(sample, elementId);
                if (overlapCount != 0) {
                    elementIdsToKeep.add(elementId);
                }

            }
        }
        assert elementIdsToKeep.size() > 0 : "kept elements cannot be null. ";
        // determine upper quartile count in each sample:
        for (final String sampleId : samplesToNormalize) {
            final DoubleArrayList countValues = new DoubleArrayList();

            for (final MutableString elementId : elementIdsToKeep) {
                countValues.add(calculator.getOverlapCount(sampleId, elementId));
            }
            Collections.sort(countValues);
            final double upperQuartile = countValues.getDouble((int) (countValues.size() * percentile));
            normalizationFactors.put(sampleId, upperQuartile);
        }
        // determine total counts over all samples considered:
        long sumOverSamples = 0;
        double sumFactors = 0;
        for (final String sampleId : samplesToNormalize) {
            sumOverSamples += calculator.getSumOverlapCounts(sampleId);
            sumFactors += normalizationFactors.get(sampleId);
        }

        // adjust the normalization factors by a constant proportion (adjustmentRatio), to bring their sum
        // to equal the sum of counts over all samples (sumOverSamples)
        final double adjustmentRatio = ((double) sumOverSamples) / sumFactors;
        for (final String sampleId : samplesToNormalize) {
            double adjustedFactor = normalizationFactors.get(sampleId);
            adjustedFactor *= adjustmentRatio;
            // force normalization factor to be at least one (to prevent divisions by zero when counts are overall
            // very small across all samples and the 75 percentile is zero.)
            final double notZero = Math.max(1, adjustedFactor);
            normalizationFactors.put(sampleId, notZero);
            LOG.info(String.format("normalization denominator %g for sample %s", notZero, sampleId));
        }

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public final double getDenominator(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                       final String sampleId) {
        return normalizationFactors.getDouble(sampleId);
    }


}
