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

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.inference.TTest;
import org.apache.commons.math.stat.inference.TTestImpl;

/**
 * Calculates T-Test P-values (two-tailed, equal variance). The T-test assess how likely the
 * mean of ln1p of the RPKMs in the first group differ from the same mean estimated in the
 * second group (requires exactly two groups with at least 2 samples in each group).
 * T-Test is applied to the ln1p, that is the natural log of the RPKM plus one.
 *
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 7:06:31 PM
 */
public class TTestCalculator extends StatisticCalculator {
    private final TTest mathCommonsTTest = new TTestImpl();

    public TTestCalculator(final DifferentialExpressionResults results) {
        this();
        setResults(results);
    }

    public TTestCalculator() {
        super();
    }

    @Override
   public boolean canDo(final String[] group) {
        return group.length == 2;
    }

    @Override
  public  DifferentialExpressionInfo evaluate(final DifferentialExpressionCalculator differentialExpressionCalculator,
                                        final NormalizationMethod method, final DifferentialExpressionResults results,
                                        final DifferentialExpressionInfo info,
                                        final String... group) {
        final String groupA = group[0];
        final String groupB = group[1];
        final ObjectArraySet<String> samplesA = differentialExpressionCalculator.getSamples(groupA);
        final ObjectArraySet<String> samplesB = differentialExpressionCalculator.getSamples(groupB);
        if (samplesA.size() < 2 || samplesB.size() < 2) {
            return info;
        }


        final int tTestStatIndex = defineStatisticId(results, "t-test", method);
        final int tStatisticStatIndex = defineStatisticId(results, "t-statistic", method);


        final double[] valuesA = new double[samplesA.size()];
        final double[] valuesB = new double[samplesB.size()];


        int i = 0;
        for (final String sample : samplesA) {
            valuesA[i++] = StrictMath.log1p(differentialExpressionCalculator.getNormalizedExpressionValue(sample, method, info.elementId));
        }

        i = 0;
        for (final String sample : samplesB) {
            valuesB[i++] = StrictMath.log1p(differentialExpressionCalculator.getNormalizedExpressionValue(sample, method, info.elementId));
        }

        double pValue = 0;
        double tStatistic = 0;

        try {
            pValue = mathCommonsTTest.homoscedasticTTest(valuesA, valuesB);
            tStatistic = mathCommonsTTest.t(valuesA, valuesB);
        } catch (MathException e) {
            pValue = Double.NaN;
        }

        info.statistics.size(results.getNumberOfStatistics());
        info.statistics.set(tTestStatIndex, pValue);
        info.statistics.set(tStatisticStatIndex, tStatistic);

        return info;
    }
}
