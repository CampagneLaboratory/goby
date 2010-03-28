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

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
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
     * Define the name of a statistic if the name was not previously defined.
     *
     * @param results     DifferentialExpressionResults instance
     * @param statisticId name of the statistic
     * @return The index of the defined statistic.
     */
    public int defineStatisticId(final DifferentialExpressionResults results, final MutableString statisticId) {
        if (!results.isStatisticDefined(statisticId)) {
            int index = results.declareStatistic(statisticId);
            statisticIds.add(statisticId);
            return index;
        }
        return results.getStatisticIndex(statisticId);
    }

    /**
     * Define the name of a statistic if the name was not previously defined. Use the abbreviation of the normalization method
     * to indicate which method was used during the calculation of this statistic.
     *
     * @param results             DifferentialExpressionResults instance
     * @param statisticId         name of the statistic
     * @param normalizationMethod Normalization method
     * @return The index of the defined statistic.
     */
    public int defineStatisticId(final DifferentialExpressionResults results, final String statisticId,
                                 final NormalizationMethod normalizationMethod) {
        return defineStatisticId(results, new MutableString(statisticId + "(" + normalizationMethod.getAbbreviation() + ")"));
    }

    /**
     * Define the name of a statistic if the name was not previously defined. Use the abbreviation of the normalization method
     * to indicate which method was used during the calculation of this statistic.
     *
     * @param results             DifferentialExpressionResults instance
     * @param statisticId         name of the statistic
     * @param normalizationMethod Normalization method
     * @return The index of the defined statistic.
     */
    public int defineStatisticId(final DifferentialExpressionResults results, final MutableString statisticId,
                                 final NormalizationMethod normalizationMethod) {
        return defineStatisticId(results, statisticId.toString(), normalizationMethod);
    }

    /**
     * Define the name of a statistic if the name was not previously defined.
     *
     * @param results     DifferentialExpressionResults instance
     * @param statisticId name of the statistic
     * @return The index of the defined statistic.
     */
    public int defineStatisticId(final DifferentialExpressionResults results, final String statisticId) {
        return defineStatisticId(results, new MutableString(statisticId));
    }

    /**
     * The name of the statistic implemented by this calculator.
     * May be null if the statistic name depends on the group name
     */
    protected ObjectList<MutableString> statisticIds=new ObjectArrayList<MutableString>();

    /**
     * Create a new StatisticCalculator with the given name.
     *
     * @param statisticIds List of statistic ids that this calculator will generate, or none if the statistic ids will
     *                     be determined by the evaluate method.
     */
    protected StatisticCalculator(final String... statisticIds) {
        for (final String id : statisticIds) {
            this.statisticIds.add(new MutableString(id));
        }

    }

    public void setResults(final DifferentialExpressionResults results) {
        this.results = results;
        for (final MutableString id : statisticIds) {
            results.declareStatistic(id);
        }
    }

    /**
     * Override this method to support a certain number of groups.
     *
     * @param group
     * @return
     */
    public abstract boolean canDo(String[] group);

    /**
     * Some implementations may not be installed on the local machine. In such cases, this method will return false and
     * client code should not call evaluate.
     * @return True when the implementation is installed on the machine.
     */
    public boolean installed(){ return true;}

    public abstract DifferentialExpressionInfo evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                                 final NormalizationMethod method,
                                                 final DifferentialExpressionResults results,
                                                 final DifferentialExpressionInfo info,
                                                 final String... group);

    /**
     * Evaluate a statistic on all the elements described in a differentialExpressionCalculator.
     *
     * @param differentialExpressionCalculator
     *               The deCalculator that keeps the data
     *               needed to evaluate the statistic
     * @param method
     * @param group  The set of groups for which the comparison is sought.  @return list of DifferentialExpressionInfo, one per element processed.
     */
    public DifferentialExpressionResults evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                           final NormalizationMethod method,
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
            evaluate(differentialExpressionCalculator, method, results, info, group);
        }
        return inputList;
    }

    /**
     * Return the statistic id that starts with prefix.
     *
     * @param prefix prefix must match the start of statistic id
     * @return the first statistic id found that starts with the prefix.
     */
    protected MutableString getMatchingStatId(final String prefix) {
        for (final MutableString id : statisticIds) {
            if (id.startsWith(prefix)) {
                return id;
            }
        }
        return null;
    }
}
