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
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import org.apache.commons.math.stat.inference.ChiSquareTest;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;
import org.apache.commons.math.MathException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.log4j.Logger;
import com.sun.tools.example.debug.gui.Environment;

/**
 * Calculates the two-tailed chi square test P-value for an observed count difference between comparison groups
 * (can assess differences between multiple groups n>=2).
 *
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 7:06:31 PM
 */
public class ChiSquareTestCalculator extends StatisticCalculator {
    private static final Logger LOG = Logger.getLogger(ChiSquareTestCalculator.class);

    public ChiSquareTestCalculator(DifferentialExpressionResults results) {
        this();
        setResults(results);


    }

    public ChiSquareTestCalculator() {
        super("chi-square-test");
    }


    public boolean canDo(String[] group) {
        return group.length >= 2;
    }


    public DifferentialExpressionInfo evaluate(DifferentialExpressionCalculator differentialExpressionCalculator,
                                               DifferentialExpressionResults results,
                                               DifferentialExpressionInfo info,
                                               String... group) {

        // expected counts in each group, assuming the counts for the DE are spread  among the groups according to sample 
        // global count proportions
        double[] expectedCounts = new double[group.length];

        //  counts observed in each group:
        long[] observedCounts = new long[group.length];
        double[] groupProportions = new double[group.length];

        int i = 0;

        double pValue = 1;
        // estimate the sumOfCountsForDE of counts over all the samples included in any group compared.
        long sumOfCountsForDE = 0;
        int numSamples = 0;
        double sumObservedCounts = 0;

        for (String oneGroupId : group) {
            ObjectArraySet<String> samplesForGroup = differentialExpressionCalculator.getSamples(oneGroupId);

            for (String sample : samplesForGroup) {

                final long observedCount = differentialExpressionCalculator.getOverlapCount(sample, info.elementId);
                final double sampleProportion = differentialExpressionCalculator.getSampleProportion(sample);
                observedCounts[i] += observedCount;
                groupProportions[i] += sampleProportion;
                sumObservedCounts += observedCount;
                numSamples++;
            }
            if (observedCounts[i] == 0) {
                // Chi Square is not defined if any observed counts are zero.
                info.statistics.size(results.getNumberOfStatistics());
                info.statistics.set(results.getStatisticIndex(STATISTIC_ID), Double.NaN);
                return info;
            }
            ++i;
        }

        i = 0;
        double nGroups = group.length;
        for (int groupIndex = 0; groupIndex < nGroups; groupIndex++) {
            expectedCounts[groupIndex] += groupProportions[groupIndex] * sumObservedCounts;
        }


        ChiSquareTest chisquare = new ChiSquareTestImpl();

        try

        {
            final double pValueRaw = chisquare.chiSquareTest(expectedCounts, observedCounts);
            // math commons can return negative p-values?
            pValue = Math.abs(pValueRaw);
        }

        catch (MaxIterationsExceededException e)

        {
            LOG.error("elementId:" + info.elementId);
            LOG.error("expected:" + DoubleArrayList.wrap(expectedCounts).toString());
            LOG.error("observed:" + LongArrayList.wrap(observedCounts).toString());
            LOG.error(e);
            pValue = 1;
        }

        catch (MathException e) {
            e.printStackTrace();

        }

        info.statistics.size(results.getNumberOfStatistics());
        info.statistics.set(results.getStatisticIndex(STATISTIC_ID), pValue);


        return info;
    }


}