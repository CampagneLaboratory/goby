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
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.math.stat.inference.TTest;
import org.apache.commons.math.stat.inference.TTestImpl;
import org.apache.commons.math.MathException;

/**
 * Calculates fold change from first group to second group (requires exactly two groups).
 *
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 7:06:31 PM
 */
public class TTestCalculator extends StatisticCalculator {

    public TTestCalculator(DifferentialExpressionResults results) {
        this();
        setResults(results);


    }

    public TTestCalculator() {
        super("t-test");
    }


    boolean canDo(String[] group) {
        return group.length == 2;
    }

    private TTest mathCommonsTTest = new TTestImpl();

    DifferentialExpressionInfo evaluate(DifferentialExpressionCalculator differentialExpressionCalculator,
                                        DifferentialExpressionResults results,
                                        DifferentialExpressionInfo info,
                                        String... group) {
        String groupA = group[0];
        String groupB = group[1];

        ObjectArraySet<String> samplesA = differentialExpressionCalculator.getSamples(groupA);
        ObjectArraySet<String> samplesB = differentialExpressionCalculator.getSamples(groupB);

        double[] valuesA = new double[samplesA.size()];
        double[] valuesB = new double[samplesB.size()];


        int i = 0;
        for (String sample : samplesA) {
            valuesA[i++] = Math.log1p(differentialExpressionCalculator.getRPKM(sample, info.elementId));
        }

        i = 0;
        for (String sample : samplesB) {
            valuesB[i++] = Math.log1p(differentialExpressionCalculator.getRPKM(sample, info.elementId));
        }

        double pValue = 0;
        double tStatistic = 0;
        final MutableString statName = new MutableString("t-statistic");
        if (!results.isStatisticDefined(statName)) {
            results.declareStatistic(statName);
        }
        try {
            pValue = mathCommonsTTest.tTest(valuesA, valuesB);
            tStatistic = mathCommonsTTest.t(valuesA, valuesB);
        } catch (MathException e) {
            pValue = Double.NaN;
        }
        info.statistics.size(results.getNumberOfStatistics());
        info.statistics.set(results.getStatisticIndex(STATISTIC_ID), pValue);
        info.statistics.set(results.getStatisticIndex(statName), tStatistic);

        return info;
    }

    final double LOG_2 = Math.log(2);

    /**
     * Calculate the log2 of x +1.
     *
     * @param x
     * @return log2(x+1)=Math.log1p(x)/Math.log(2)
     */
    private double log2(final double x) {

        return Math.log1p(x) / LOG_2;
    }

}