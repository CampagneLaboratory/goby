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

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;

import java.io.PrintWriter;
import java.util.List;

/**
 * Store a list of results from differential expression statistics.
 *
 * @author Fabien Campagne
 *         Date: Jan 12, 2010
 *         Time: 5:34:36 PM
 */
public class DifferentialExpressionResults extends ObjectArrayList<DifferentialExpressionInfo> {
    private final IndexedIdentifier statisticIds = new IndexedIdentifier();
    private final List<MutableString> sortedStatisticIds = new ObjectArrayList<MutableString>();
    private IntArrayList averageCountPerGroupIndexes = new IntArrayList();
    private boolean omitNonInformativeColumns = false;

    /**
     * Declare a new statistic.
     *
     * @param statisticId
     */
    public synchronized void declareStatistic(final String statisticId) {
        declareStatistic(new MutableString(statisticId));
    }

    /**
     * Declare a new statistic.
     *
     * @param statisticId Identifier for the new statistic.
     * @return the index of the statistic.
     */
    public int declareStatistic(final MutableString statisticId) {
       int index= statisticIds.registerIdentifier(statisticId);
       // TODO There is no reason for this class to know anything about the name of columns it may contain. Remove this code.
        if (statisticId.startsWith("average count group ")) {
            averageCountPerGroupIndexes.add(sortedStatisticIds.size());
        }
        sortedStatisticIds.add(statisticId);
        return index;
    }

    public boolean isOmitNonInformativeColumns() {
        return omitNonInformativeColumns;
    }

    public void setOmitNonInformativeColumns(boolean omitNonInformativeColumns) {
        this.omitNonInformativeColumns = omitNonInformativeColumns;
    }

    public double getStatistic(final DifferentialExpressionInfo info, final MutableString statisticId) {
        return info.statistics.get(statisticIds.get(statisticId));
    }

    public int getStatisticIndex(final MutableString statisticId) {
        return statisticIds.getInt(statisticId);
    }

    public int getNumberOfStatistics() {
        return statisticIds.size();
    }

    @Override
    public String toString() {
        final MutableString buffer = new MutableString();
        buffer.append("element-id ");
        for (final MutableString statId : sortedStatisticIds) {
            buffer.append(statId);
            buffer.append(" ");
        }
        buffer.append("\n");
        for (final DifferentialExpressionInfo info : this) {
            buffer.append(info.toString());
            buffer.append('\n');
        }

        return buffer.toString();
    }

    public boolean isStatisticDefined(final MutableString groupId) {
        return statisticIds.containsKey(groupId);
    }

    public synchronized IntArrayList getAverageCountPerGroupIndexes() {
        return averageCountPerGroupIndexes;
    }

    /**
     * Write results with delimiter.
     *
     * @param printWriter
     * @param delimiter
     */
    public void write(final PrintWriter printWriter, final char delimiter) {
        InformativeColumns informativeColumns = null;
        if (omitNonInformativeColumns) {
            informativeColumns = new InformativeColumns(sortedStatisticIds.size(), new InformativeNonZeroNonNaN());
            for (final DifferentialExpressionInfo info : this.subList(0, size())) {
                if (info.checkInformativeColumns(informativeColumns)) {
                    break;
                }
            }
        }

        printWriter.append("element-id");
        for (int i = 0; i < sortedStatisticIds.size(); i++) {
            if (informativeColumns == null || informativeColumns.isColumnInformative(i)) {
                printWriter.append(delimiter);
                printWriter.append(sortedStatisticIds.get(i));
            }
        }
        printWriter.append("\n");

        for (final DifferentialExpressionInfo info : this.subList(0, size())) {
            if (info.informative(getAverageCountPerGroupIndexes())) {
                info.write(printWriter, delimiter, informativeColumns);
                printWriter.append("\n");
            }
        }
        printWriter.flush();
    }

    public int getStatisticIndex(final String statisticId) {
        return getStatisticIndex(new MutableString(statisticId));
    }
}
