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

import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;

/**
 * Calculate differential expression statistics for lists of elements under study
 * (e.g., genes, exons).
 *
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 7:02:02 PM
 */
public abstract class StatisticCalculator {
    protected DifferentialExpressionResults results;

    /**
     * The name of the statistic implemented by this calculator.
     * May be null if the statistic name depends on the group name
     */
    protected final MutableString statisticId;

    /**
     * Create a new StatisticCalculator with the given name.
     * @param statisticId
     */
    protected StatisticCalculator(final String statisticId) {
        if (statisticId != null) {
            this.statisticId = new MutableString(statisticId);
        } else {
            this.statisticId = null;
        }
    }

    public void setResults(final DifferentialExpressionResults results) {
        this.results = results;
        if (statisticId != null) {
            results.declareStatistic(statisticId);
        }
    }

    /**
     * Override this method to support a certain number of groups.
     *
     * @param group
     * @return
     */
    abstract boolean canDo(String[] group);

    abstract DifferentialExpressionInfo evaluate(DifferentialExpressionCalculator differentialExpressionCalculator,
                                                 DifferentialExpressionResults results,
                                                 DifferentialExpressionInfo info,
                                                 String... group);

    /**
     * Evaluate a statistic on all the elements described in a differentialExpressionCalculator.
     *
     * @param differentialExpressionCalculator The deCalculator that keeps the data
     * needed to evaluate the statistic
     * @param group The set of groups for which the comparison is sought.
     * @return list of DifferentialExpressionInfo, one per element processed.
     */
    DifferentialExpressionResults evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                           final DifferentialExpressionResults inputList,
                                           final String... group) {
        if (inputList.size() == 0) {
            // define an info elements for each element ids defined by the calculator:
            final ObjectSet<MutableString> elementIds = differentialExpressionCalculator.getElementIds();
            for (final MutableString elementId : elementIds) {
                final DifferentialExpressionInfo info = new DifferentialExpressionInfo(elementId);
                results.add(info);
            }
        }

        for (final DifferentialExpressionInfo info : inputList) {
            evaluate(differentialExpressionCalculator, results, info, group);
        }
        return inputList;
    }
}
