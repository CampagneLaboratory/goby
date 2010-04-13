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

/**
 * Assists with determining which columns in a large grid of "double" data are informative.
 * Optimized to not double-check columns that are known to be informative and when all columns
 * are known to be informative it does no more checking at all.
 *
 * Informative is here defined as not NaN and > 0.
 *
 * @author Kevin Dorff
 */
public class InformativeColumns {

    /**
     * True if all columns have been determined to be informative. This will only be
     * set after the last entry of a row has been checked.
     */
    public boolean allColumnsInformative;

    /**
     * The number of columns in each row of the dataset.
     */
    int numberOfColumns;
    /**
     * The current position in the current row we are checking. Automatically returns to
     * 0 when we have checked numberOfColumns values.
     */
    int positionInCurrentRow;

    /**
     * The data, one boolean per numberOfColumns. Starts off as all false since rows are not
     * yet known to be informative.
     */
    boolean[] data;

    /**
     * The object that defines what is informative.
     */
    public InformativeDouble informativeDouble;

    /**
     * Create an InformativeColumns object for numberOfColumns.
     * @param numberOfColumns the fixed number of columns in the dataset for any number of rows.
     * @param informativeDouble the object that defines what is informative.
     */
    public InformativeColumns(final int numberOfColumns, final InformativeDouble informativeDouble) {
        assert numberOfColumns > 0;
        this.informativeDouble = informativeDouble;
        this.numberOfColumns = numberOfColumns;
        allColumnsInformative = false;
        positionInCurrentRow = 0;
        data = new boolean[numberOfColumns];
    }

    /**
     * An entire row of data has been checked. Check if all columns are now informative.
     */
    private void newRowCheckAllInformative() {
        if (allColumnsInformative) {
            return;
        }
        allColumnsInformative = true;
        for (final boolean column : data) {
            if (!column) {
                allColumnsInformative = false;
            }
        }
        positionInCurrentRow = 0;
    }

    /**
     * Return true if all columns are informative.
     * @return true if all columns are informative.
     */
    public boolean isAllColumnsInformative() {
        return allColumnsInformative;
    }

    /**
     * Check if the next value is informative.
     * @param value the value to check
     */
    public void checkInformative(final double value) {
        if (allColumnsInformative) {
            return;
        }
        if (!data[positionInCurrentRow]) {
            if (informativeDouble.isInformative(value)) {
                data[positionInCurrentRow] = true;
            }
        }
        if (++positionInCurrentRow == numberOfColumns) {
            newRowCheckAllInformative();
        }
    }

    /**
     * Check if the specified column is informative.
     * @param column the column to check
     * @return true if that column is informative
     */
    public boolean isColumnInformative(final int column) {
        return data[column];
    }
}
