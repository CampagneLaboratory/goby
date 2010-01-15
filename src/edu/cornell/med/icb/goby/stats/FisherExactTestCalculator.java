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
 * Calculates Fisher exact test P-value for an observed count difference between comparison groups (requires exactly two groups).
 *
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 7:06:31 PM
 */
public class FisherExactTestCalculator extends StatisticCalculator {

    public FisherExactTestCalculator(DifferentialExpressionResults results) {
        this();
        setResults(results);


    }

    public FisherExactTestCalculator() {
        super("fisher-exact-test");
    }


    boolean canDo(String[] group) {
        return group.length == 2;
    }


    DifferentialExpressionInfo evaluate(DifferentialExpressionCalculator differentialExpressionCalculator,
                                        DifferentialExpressionResults results,
                                        DifferentialExpressionInfo info,
                                        String... group) {
        String groupA = group[0];
        String groupB = group[1];

        ObjectArraySet<String> samplesA = differentialExpressionCalculator.getSamples(groupA);
        ObjectArraySet<String> samplesB = differentialExpressionCalculator.getSamples(groupB);

        int sumCountInA = 0;// = new double[samplesA.size()];
        int sumCountInB = 0;// = new double[samplesB.size()];


        for (String sample : samplesA) {
            sumCountInA += differentialExpressionCalculator.getOverlapCount(sample, info.elementId);
        }


        for (String sample : samplesB) {
            sumCountInB += differentialExpressionCalculator.getOverlapCount(sample, info.elementId);
        }
        int totalCountInA = 0;
        int totalCountInB = 0;


        for (String sample : samplesA) {
            totalCountInA += differentialExpressionCalculator.getSumOverlapCounts(sample);
        }
        for (String sample : samplesB) {
            totalCountInB += differentialExpressionCalculator.getSumOverlapCounts(sample);
        }


        double pValue = 1;

        /**
         * Calculates the 2-tailed gominer.Fisher p value, using the counts that are
         * expected to be readily available in the system
         *
         * @param totalChanged The number of genes changed in the experiment
         * N1
         * @param changedInNode The number of genes changed in a particular node
         * x
         * @param total The total number of genes in the experiment (changed and
         * unchanged)     N1+N2
         * @param inNode The number of genes in the node (changed and unchanged)
         * x+y
         * @return 2-tailed gominer.Fisher p value
         */
        // public double fisher(final int totalChanged, final int changedInNode, final int total, final int inNode) {
        Fisher fisher = new Fisher();

        pValue = fisher.fisher(totalCountInA, sumCountInA, totalCountInA + totalCountInB, sumCountInA + sumCountInB);

       
        /* Test : fisher.fisher(40,10,100,30)=
                     Fisher's Exact Test
        http://www.langsrud.com/fisher.htm
        ------------------------------------------
         TABLE = [ 10 , 20 , 30 , 40 ]
        Left   : p-value = 0.2533310713617698
        Right  : p-value = 0.8676419647894328
        2-Tail : p-value = 0.5044757698516504
        ------------------------------------------
        */
        info.statistics.size(results.getNumberOfStatistics());
        info.statistics.set(results.getStatisticIndex(STATISTIC_ID), pValue);


        return info;
    }


}