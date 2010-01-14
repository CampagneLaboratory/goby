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

/**
 * Calculates fold change from first group to second group (requires exactly two groups).
 *
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 7:06:31 PM
 */
public class FoldChangeCalculator extends StatisticCalculator {

    public FoldChangeCalculator(DifferentialExpressionResults results) {
        this();
        setResults(results);


    }

    public FoldChangeCalculator() {
        super("fold-change");
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
        double averageA = 0;
        double averageB = 0;


        for (String sample : samplesA) {
            averageA += differentialExpressionCalculator.getRPKM(sample, info.elementId);
        }
        for (String sample : samplesB) {
            averageB += differentialExpressionCalculator.getRPKM(sample, info.elementId);
        }

        double foldChangeStatistic = averageA / averageB;
        info.statistics.size(results.getNumberOfStatistics());
        info.statistics.set(results.getStatisticIndex(STATISTIC_ID), foldChangeStatistic);

        return info;
    }


}
