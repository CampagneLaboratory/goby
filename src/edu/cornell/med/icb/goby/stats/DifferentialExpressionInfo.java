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

import java.io.PrintWriter;

/**
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 6:57:49 PM
 */
public class DifferentialExpressionInfo {
    final MutableString elementId;
    final DoubleArrayList statistics = new DoubleArrayList();

    public DifferentialExpressionInfo(final String elementId) {
        super();
        this.elementId = new MutableString(elementId);
    }

    public DifferentialExpressionInfo(final MutableString elementId) {
        super();
        this.elementId = elementId;
    }

    @Override
    public String toString() {
        return String.format("[ %s %s ]", elementId, statistics.toString());
    }

    public void write(final PrintWriter printWriter, final char delimiter, final InformativeColumns informativeColumns) {
        printWriter.append(elementId);
        for (int i = 0; i < statistics.size(); i++) {
            if (informativeColumns == null || informativeColumns.isColumnInformative(i)) {
                final double value = statistics.get(i);
                printWriter.append(delimiter);
                printWriter.append(String.format("%g", value));
            }
        }
    }

    public void checkInformativeColumns(final InformativeColumns informativeColumns) {
        if (informativeColumns.isAllColumnsInformative()) {
            return;
        }
        for (double value : statistics) {
            informativeColumns.checkInformative(value);
        }
    }

    /**
     * Is the DE informative?
     * @param averageCountPerGroupIndexes list of average counts for group. If this isn't null,
     * the value for at least one of these indexes must be informative (!NaN and > 0).
     * @return if this line of data is informative
     */
    public boolean informative(final IntArrayList averageCountPerGroupIndexes) {
        if (averageCountPerGroupIndexes != null) {
            boolean atLeastOneGroupAverageNotZero = false;
            for (int informativeRequiredIndex : averageCountPerGroupIndexes) {
                double value = statistics.get(informativeRequiredIndex);
                if (!Double.isNaN(value) && value > 0) {
                    atLeastOneGroupAverageNotZero = true;
                }
            }
            if (!atLeastOneGroupAverageNotZero) {
                return false;
            }
        }

        boolean informative = false;
        for (final double value : statistics) {
            if (!Double.isNaN(value) && value > 0) {
                // require something else than NaN or zero to be have an informative DE.
                informative = true;
            }
        }
        return informative;
    }
}
