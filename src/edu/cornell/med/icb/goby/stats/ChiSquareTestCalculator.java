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
import org.apache.commons.math.stat.inference.ChiSquareTest;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;
import org.apache.commons.math.MathException;
import gominer.Fisher;

/**
 * Calculates the two-tailed chi square test P-value for an observed count difference between comparison groups
 * (can assess differences between multiple groups n>=2).
 *
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 7:06:31 PM
 */
public class ChiSquareTestCalculator extends StatisticCalculator {

    public ChiSquareTestCalculator(DifferentialExpressionResults results) {
        this();
        setResults(results);


    }

    public ChiSquareTestCalculator() {
        super("chi-square-test");
    }


   public  boolean canDo(String[] group) {
        return group.length >= 2;
    }


    public DifferentialExpressionInfo evaluate(DifferentialExpressionCalculator differentialExpressionCalculator,
                                        DifferentialExpressionResults results,
                                        DifferentialExpressionInfo info,
                                        String... group) {

        //    totalCountInA, totalCountInB,...
        double[] expectedCounts = new double[group.length];


        //  sumCountInA, sumCountInB ,...
        long[] observedCounts = new long[group.length];

        int i = 0;
        for (String groupI : group) {
            ObjectArraySet<String> samplesI = differentialExpressionCalculator.getSamples(groupI);

            for (String sample : samplesI) {
                observedCounts[i] += differentialExpressionCalculator.getOverlapCount(sample, info.elementId);
            }

            for (String sample : samplesI) {
                expectedCounts[i] += differentialExpressionCalculator.getSumOverlapCounts(sample);
            }
            ++i;

        }
        double pValue;

            ChiSquareTest chisquare = new ChiSquareTestImpl();

            try {
                pValue = chisquare.chiSquareTest(expectedCounts, observedCounts);
            } catch (MathException e) {
                pValue = Double.NaN;
            }
       
        info.statistics.size(results.getNumberOfStatistics());
        info.statistics.set(results.getStatisticIndex(STATISTIC_ID), pValue);


        return info;
    }


}