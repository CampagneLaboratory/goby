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

import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntCollection;
import it.unimi.dsi.lang.MutableString;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Formatter;

/**
 * Helper class to write tab delimited statistic result files.
 *
 * @author Fabien Campagne
 *         Date: Sep 23, 2010
 *         Time: 1:26:56 PM
 */
public class TSVWriter {
    PrintWriter outWriter;

    public TSVWriter(PrintWriter outWriter) {
        this.outWriter = outWriter;
        columnIds = new IndexedIdentifier();
        columnIds.defaultReturnValue(COLUMN_NOT_DEFINED);

    }

    IndexedIdentifier columnIds = new IndexedIdentifier();
    Formatter formatter;

    /**
     * Define a new column, built with String.format and optional parameters.
     *
     * @param columnId
     * @param option
     * @return
     */
    public int defineColumn(String columnId, String... option) {
        String id = String.format(columnId, option);
        final int index = columnIds.registerIdentifier(new MutableString(id));
        values = new CharSequence[columnIds.size()];
        return index;
    }

    /**
     * Define new columns, built with String.format and optional parameters.
     * Columns are defined in groups of ids with different options for the same id following each other.
     *
     * @param options set of String.format ids parameters (exactly one per id)
     * @param ids     Set of identifier format to create new column identifiers.
     */
    public void defineColumnSet(String[] options, String... ids) {

        for (String columnId : ids) {   //define groups of columns, grouped by option format
            for (String option : options) {
                String id = String.format(columnId, option);
                columnIds.registerIdentifier(new MutableString(id));
                values = new CharSequence[columnIds.size()];
            }
        }


    }

    public static int COLUMN_NOT_DEFINED = -1;

    public void close() {
        outWriter.close();
    }

    CharSequence[] values;

    public void setValue(int columnIndex, CharSequence value) {
        if (columnIndex == COLUMN_NOT_DEFINED) return;
        values[columnIndex] = value;
    }

    public void setValue(String value, String id, String... option) {
        MutableString idMut = new MutableString(String.format(id, option));
        int columnIndex = this.columnIds.getInt(idMut);
        if (columnIndex == COLUMN_NOT_DEFINED) {
            System.out.println("col not found: " + idMut);
            return;
        }
        values[columnIndex] = value;
    }

    public void setValue(int value, String id, String... option) {
        setValue(String.format("%d", value), id, option);
    }

    public void setValue(double value, String id, String... option) {
        setValue(String.format("%g", value), id, option);
    }

    public void setValue(int columnIndex, double value) {
        if (columnIndex == COLUMN_NOT_DEFINED) return;
        values[columnIndex] = String.format("%g", value);
    }

    public void setValue(int columnIndex, int value) {
        if (columnIndex == COLUMN_NOT_DEFINED) return;
        values[columnIndex] = String.format("%d", value);
    }

    public void setValue(int columnIndex, char value) {
        if (columnIndex == COLUMN_NOT_DEFINED) return;
        values[columnIndex] = String.format("%c", value);
    }

    public void setValue(int columnIndex, byte value) {
        if (columnIndex == COLUMN_NOT_DEFINED) return;
        values[columnIndex] = String.format("%d", value);
    }

    public void setValue(int columnIndex, long value) {
        if (columnIndex == COLUMN_NOT_DEFINED) return;
        values[columnIndex] = String.format("%d", value);
    }

    public void writeRecord() {
        int index = 0;
        int lastIndex = values.length - 1;
        for (CharSequence value : values) {
            outWriter.print(value);
            if (index != lastIndex) {
                outWriter.print("\t");
            } else {
                outWriter.println();
            }

            ++index;

        }
        Arrays.fill(values, "");
    }

    public void writeHeader() {
        final IntCollection intCollection = columnIds.values();

        int max = -1;
        for (int columnIndex : intCollection) {
            max = Math.max(columnIndex, max);
        }
        int index = 0;
        int lastIndex = 0;
        DoubleIndexedIdentifier reverse = new DoubleIndexedIdentifier(columnIds);
        for (int columnIndex = 0; columnIndex <= max; columnIndex++) {
            outWriter.print(reverse.getId(columnIndex));
            if (index != max) {
                outWriter.print("\t");
            } else {
                outWriter.println();
            }

            ++index;

        }

    }

    /**
     * This method does nothing in TSVWriter.
     *
     * @param shortId
     * @param numValues
     * @param type
     * @param description
     * @param columnId
     * @param option
     */
    public void defineColumnAttributes(String shortId, int numValues, ColumnType type, String description, String columnId, String... option) {

    }

    /**
     * This method does nothing in TSVWriter.
     */
    public void defineColumnAttributes(int numValues, ColumnType type, String[] options, String... columnIds) {
       
    }
}
