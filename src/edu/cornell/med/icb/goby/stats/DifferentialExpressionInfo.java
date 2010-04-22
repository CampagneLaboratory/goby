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

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.lang.MutableString;

import java.io.IOException;
import java.io.Writer;

/**
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 6:57:49 PM
 */
public class DifferentialExpressionInfo {
    private final MutableString elementId;
    final DoubleArrayList statistics = new DoubleArrayList();
    private final InformativeDouble informativeDouble = new InformativeNonZeroNonNaN();

    public DifferentialExpressionInfo(final String elementId) {
        super();
        this.elementId = new MutableString(elementId);
    }

    public DifferentialExpressionInfo(final MutableString elementId) {
        super();
        this.elementId = elementId;
    }

    public MutableString getElementId() {
        return elementId;
    }

    @Override
    public String toString() {
        return String.format("[ %s %s ]", elementId, statistics.toString());
    }

    public void write(final Writer writer, final char delimiter,
                      final InformativeColumns informativeColumns,
                      final DifferentialExpressionCalculator deCalculator) throws IOException {
        writer.append(elementId);
        writer.append(delimiter);
        final String elementType = deCalculator.getElementType(elementId).toString();

        writer.append(elementType);
        for (int i = 0; i < statistics.size(); i++) {
            if (informativeColumns == null || informativeColumns.isColumnInformative(i)) {
                writer.append(delimiter);
                final double value = statistics.get(i);
                writer.append(String.format("%g", value));
            }
        }
    }

    /**
     * Check the data in the current row to see which colums are informative.
     * @param informativeColumns the object that helps track which columns are informative
     * @return true if all columns have been determined to be informative
     */
    public boolean checkInformativeColumns(final InformativeColumns informativeColumns) {
        for (final double value : statistics) {
            informativeColumns.checkInformative(value);
        }
        return informativeColumns.isAllColumnsInformative();
    }

    /**
     * Is the DE informative?
     * @return if this line of data is informative
     */
    public boolean informative() {
        return informative(null);
    }

    /**
     * Is the DE informative?
     * @param averageCountPerGroupIndexes list of average counts for group. If this isn't null,
     * the value for at least one of these indexes must be informative (!NaN and > 0).
     * @return if this line of data is informative
     */
    public boolean informative(final IntArrayList averageCountPerGroupIndexes) {
        if (averageCountPerGroupIndexes != null && averageCountPerGroupIndexes.size() > 0) {
            boolean atLeastOneGroupAverageNotZero = false;
            for (final int informativeRequiredIndex : averageCountPerGroupIndexes) {
                final double value = statistics.get(informativeRequiredIndex);
                if (informativeDouble.isInformative(value)) {
                    atLeastOneGroupAverageNotZero = true;
                }
            }
            if (!atLeastOneGroupAverageNotZero) {
                return false;
            }
        }

        boolean informative = false;
        for (final double value : statistics) {
            if (informativeDouble.isInformative(value)) {
                // require something else than NaN or zero to be have an informative DE.
                informative = true;
            }
        }
        return informative;
    }
}
