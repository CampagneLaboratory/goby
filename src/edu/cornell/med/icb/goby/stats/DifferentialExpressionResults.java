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
import it.unimi.dsi.lang.MutableString;
import edu.cornell.med.icb.identifier.IndexedIdentifier;

import java.io.PrintWriter;

/**
 * Store a list of results from differential expression statistics.
 * 
 * @author Fabien Campagne
 *         Date: Jan 12, 2010
 *         Time: 5:34:36 PM
 */
public class DifferentialExpressionResults extends ObjectArrayList<DifferentialExpressionInfo> {

    private IndexedIdentifier statisticIds = new IndexedIdentifier();
    private ObjectArrayList<MutableString> sortedStatisticIds = new ObjectArrayList<MutableString>();

    /**
     * Declare a new statistic.
     *
     * @param statisticId
     */
    public void declareStatistic(String statisticId) {
        final MutableString mutableString = new MutableString(statisticId);
        statisticIds.registerIdentifier(mutableString);
        sortedStatisticIds.add(mutableString);
    }

    /**
     * Declare a new statistic.
     *
     * @param statisticId
     */
    public void declareStatistic(MutableString statisticId) {
        statisticIds.registerIdentifier(statisticId);
        sortedStatisticIds.add(statisticId);
    }

    public double getStatistic(DifferentialExpressionInfo info, MutableString statisticId) {
        return info.statistics.get(statisticIds.get(statisticId));
    }

    public int getStatisticIndex(MutableString statisticId) {
        return statisticIds.getInt(statisticId);
    }

    public int getNumberOfStatistics() {
        return statisticIds.size();
    }

    @Override
    public String toString() {

        MutableString buffer = new MutableString();
        buffer.append("element-id ");
        for (MutableString statId : sortedStatisticIds) {
            buffer.append(statId);
            buffer.append(" ");
        }
        buffer.append("\n");
       for (DifferentialExpressionInfo info: this) {
            buffer.append(info.toString());
            buffer.append('\n'); 
       }

        return buffer.toString();
    }

    public boolean isStatisticDefined(MutableString groupId) {
        return statisticIds.containsKey(groupId);
    }

    /**
     * Write results with delimiter.
     *
     * @param printWriter
     * @param delimiter
     */
    public void write(PrintWriter printWriter, char delimiter) {

        printWriter.append("element-id ");
        for (MutableString statId : sortedStatisticIds) {
            printWriter.append(delimiter);
            printWriter.append(statId);

        }
        printWriter.append("\n");
        for (DifferentialExpressionInfo info : this.subList(0, size())) {
            if (info.informative()) {
            info.write(printWriter, delimiter);
            printWriter.append("\n");
            }
        }
        printWriter.flush();
    }

    public int getStatisticIndex(String statisticId) {
        return getStatisticIndex(new MutableString(statisticId));
    }
}
