/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;

import java.io.IOException;
import java.io.Writer;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
    private final IntArrayList averageCountPerGroupIndexes = new IntArrayList();
    private boolean omitNonInformativeColumns;

    public DifferentialExpressionResults(int capacity) {
        super(capacity);
    }

    public DifferentialExpressionResults() {
        super();
    }

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
        final int index = statisticIds.registerIdentifier(statisticId);

        sortedStatisticIds.add(statisticId);
        return index;
    }

    /**
     * Indicate that the statistic at this index should be considered to determine informative rows.
     * @param statIndex index of the statistic to consider.
     */
    public void recordInformativeColumnIndex(int statIndex) {

       averageCountPerGroupIndexes.add(statIndex);
    }

    public boolean isOmitNonInformativeColumns() {
        return omitNonInformativeColumns;
    }

    public void setOmitNonInformativeColumns(final boolean omitNonInformativeColumns) {
        this.omitNonInformativeColumns = omitNonInformativeColumns;
    }

    public double getStatistic(final DifferentialExpressionInfo info, final MutableString statisticId) {
        return info.statistics.get(statisticIds.get(statisticId));
    }

    public double getStatistic(final DifferentialExpressionInfo info, final int statIndex) {
        return info.statistics.get(statIndex);
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

    public MutableString getStatisticIdForIndex(final int statisticIndex) {
        return sortedStatisticIds.get(statisticIndex);
    }

    public IntList statisticsIndexesFor(final String statPrefix, final NormalizationMethod method) {
        final IntList results = new IntArrayList();
        final String methodAbbr = "(" + method.getAbbreviation() + ")";
        final Set<Map.Entry<MutableString, Integer>> entrySet = statisticIds.entrySet();
        for (final Map.Entry<MutableString, Integer> entry : entrySet) {
            // First try with the method abbreviation
            if (entry.getKey().startsWith(statPrefix) && entry.getKey().endsWith(methodAbbr)) {
                results.add(entry.getValue());
            }
        }
        if (results.size() == 0) {
            // With method abbreviation matched nothing. Try without the method abbreviation.
            for (final Map.Entry<MutableString, Integer> entry : entrySet) {
                if (entry.getKey().startsWith(statPrefix)) {
                    results.add(entry.getValue());
                }
            }
        }
        return results;
    }

    public synchronized IntArrayList getAverageCountPerGroupIndexes() {
        return averageCountPerGroupIndexes;
    }

    /**
     * Write results with delimiter.
     *
     * @param writer
     * @param delimiter
     * @param deCalculator
     */
    public void write(final Writer writer, final char delimiter,
                      final DifferentialExpressionCalculator deCalculator) throws IOException {
        InformativeColumns informativeColumns = null;
        if (omitNonInformativeColumns) {
            informativeColumns = new InformativeColumns(sortedStatisticIds.size(), new InformativeNonZeroNonNaN());
            for (final DifferentialExpressionInfo info : this.subList(0, size())) {
                if (info.checkInformativeColumns(informativeColumns)) {
                    break;
                }
            }
        }

        writer.append("element-id");
        writer.append(delimiter);
        writer.append("element-type");
        for (int i = 0; i < sortedStatisticIds.size(); i++) {
            if (informativeColumns == null || informativeColumns.isColumnInformative(i)) {
                writer.append(delimiter);
                writer.append(sortedStatisticIds.get(i));
            }
        }
        writer.append("\n");

        for (final DifferentialExpressionInfo info : this.subList(0, size())) {
            if (info.informative(getAverageCountPerGroupIndexes())) {
                info.write(writer, delimiter, informativeColumns, deCalculator);
                writer.append("\n");
            }
        }
        writer.flush();
    }

    public int getStatisticIndex(final String statisticId) {
        return getStatisticIndex(new MutableString(statisticId));
    }
}
