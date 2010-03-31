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
    }

    @Override
    public boolean canDo(final String[] group) {
        return group.length == 2;
    }

    @Override
    public boolean installed() {
        final Rengine rEngine = GobyRengine.getInstance().getRengine();
        return rEngine != null && rEngine.isAlive();

    }

    @Override
    public DifferentialExpressionInfo evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                               final NormalizationMethod method, final DifferentialExpressionResults results,
                                               final DifferentialExpressionInfo info,
                                               final String... group) {
        // we can only perform the evaluation if R is running and alive.
        final Rengine rengine = GobyRengine.getInstance().getRengine();
        if (rengine != null && rengine.isAlive()) {
            final String groupA = group[0];
            final String groupB = group[1];
            // TODO correct sumCountIn? with normalization method.
            final int statIndex = defineStatisticId(results, "fisher-exact-R");

            final ObjectArraySet<String> samplesA = differentialExpressionCalculator.getSamples(groupA);
            final ObjectArraySet<String> samplesB = differentialExpressionCalculator.getSamples(groupB);

            int sumCountInA = 0;
            int sumCountInB = 0;
            // TODO correct sumCountIn? with normalization method.
            for (final String sample : samplesA) {
                sumCountInA += differentialExpressionCalculator.getOverlapCount(sample, info.elementId);
            }
            // TODO correct sumCountIn? with normalization method.
            for (final String sample : samplesB) {
                sumCountInB += differentialExpressionCalculator.getOverlapCount(sample, info.elementId);
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
        return info;
    }
}
