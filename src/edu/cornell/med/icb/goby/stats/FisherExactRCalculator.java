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

import edu.cornell.med.icb.goby.R.FisherExact;
import edu.cornell.med.icb.goby.R.GobyRengine;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.rosuda.JRI.Rengine;

/**
 * Calculates Fisher exact test P-value for an observed count difference between comparison
 * groups (requires exactly two groups).
 *
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 7:06:31 PM
 */
public class FisherExactRCalculator extends StatisticCalculator {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(FisherExactRCalculator.class);

    public FisherExactRCalculator(final DifferentialExpressionResults results) {
        this();
        setResults(results);
    }

    public FisherExactRCalculator() {
        super();
        rEngine = GobyRengine.getInstance().getRengine();
        installed = rEngine != null && rEngine.isAlive();
    }

    @Override
    public boolean canDo(final String[] group) {
        return group.length == 2;
    }

    boolean installed;

    @Override
    public boolean installed() {

        // we can only perform the evaluation if R is running and alive.

        return installed;
    }

    Rengine rEngine;

    @Override
    public DifferentialExpressionInfo evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                               final NormalizationMethod method, final DifferentialExpressionResults results,
                                               final DifferentialExpressionInfo info,
                                               final String... group) {
        synchronized (this) {
            if (installed) {
                final String groupA = group[0];
                final String groupB = group[1];

                // TODO correct sumCountIn? with normalization method.
                final int statIndex = defineStatisticId(results, "fisher-exact-R", method, group);

                final ObjectArraySet<String> samplesA = differentialExpressionCalculator.getSamples(groupA);
                final ObjectArraySet<String> samplesB = differentialExpressionCalculator.getSamples(groupB);

                int sumCountInA = 0;
                int sumCountInB = 0;
                // TODO correct sumCountIn? with normalization method.
                for (final String sample : samplesA) {
                    sumCountInA += differentialExpressionCalculator.getOverlapCount(sample, info.getElementId());
                }
                // TODO correct sumCountIn? with normalization method.
                for (final String sample : samplesB) {
                    sumCountInB += differentialExpressionCalculator.getOverlapCount(sample, info.getElementId());
                }
                int totalCountInA = 0;
                int totalCountInB = 0;


                for (final String sample : samplesA) {
                    totalCountInA += differentialExpressionCalculator.getSumOverlapCounts(sample);
                }
                for (final String sample : samplesB) {
                    totalCountInB += differentialExpressionCalculator.getSumOverlapCounts(sample);
                }

                final int sumCountNotInA = totalCountInA - sumCountInA;
                final int sumCountNotInB = totalCountInB - sumCountInB;

                // Build a contingency matrix as follows:
                //                  groupA            groupB
                // hasCounts    sumCountInA        sumCountInB
                // noCounts     sumCountNotInA     sumCountNotInB
                final FisherExact.Result result =
                        FisherExact.fexact(sumCountInA, sumCountNotInA, sumCountInB, sumCountNotInB);
                if (LOG.isDebugEnabled()) {
                    LOG.debug(result);
                }
                final double pValue = result.getPValue();
                info.statistics.size(results.getNumberOfStatistics());
                info.statistics.set(statIndex, pValue);
            }
        }

        return info;
    }

    /**
     * Estimate the Fisher P-value given the contingency table:
     * //               group0       group1
     * // condition0   count00      count01
     * // condition1   count10      count11
     *
     * @param count00
     * @param count10
     * @param count01
     * @param count11
     * @return P-value of observing a contingency table that extreme by random distribution among the cells.
     */
    public static double getFisherPValue(int count00, int count10, int count01, int count11) {
        final FisherExact.Result result =
                FisherExact.fexact(count00, count10, count01, count11);
        if (LOG.isDebugEnabled()) {
            LOG.debug(result);
        }
        final double pValue = result.getPValue();
        return pValue;
    }

    /**
     * Estimate the Fisher one-tailed lesser P-value given the contingency table:
     * //               group0       group1
     * // condition0   count00      count01
     * // condition1   count10      count11
     *
     * @param count00
     * @param count10
     * @param count01
     * @param count11
     * @return P-value of observing a contingency table that extreme by random distribution among the cells.
     */
    public static double getFisherOneTailedLesserPValue(int count00, int count10, int count01, int count11) {
        final FisherExact.Result result =
                FisherExact.fexactLesser(count00, count10, count01, count11);
        if (LOG.isDebugEnabled()) {
            LOG.debug(result);
        }
        final double pValue = result.getPValue();
        return pValue;

    }
}
