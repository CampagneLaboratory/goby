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
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.modes.GobyDriver;
import edu.cornell.med.icb.util.VersionUtils;
import edu.rit.mp.buf.SharedBooleanArrayBuf;

import java.io.PrintWriter;
import java.util.Formatter;
import java.util.Arrays;
import java.util.Collections;

import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.ints.*;

/**
 * Helper class to write tab delimited statistic result files.
 *
 * @author Fabien Campagne
 *         Date: Sep 23, 2010
 *         Time: 1:26:56 PM
 */
public class StatisticsWriter {
    PrintWriter outWriter;
    private Int2ObjectMap<ColumnType> indexTypes;
    private boolean VCFmode;
    private Int2ObjectMap<String> indexDescriptions;
    private Int2IntMap indexNumValues;
    private Int2ObjectMap<String> indexShortIds;

    public StatisticsWriter(PrintWriter outWriter) {
        this.outWriter = outWriter;
        columnIds = new IndexedIdentifier();
        columnIds.defaultReturnValue(COLUMN_NOT_DEFINED);
        indexTypes = new Int2ObjectOpenHashMap<ColumnType>();
        indexDescriptions = new Int2ObjectOpenHashMap<String>();
        indexShortIds = new Int2ObjectOpenHashMap<String>();
        indexNumValues = new Int2IntArrayMap();
    }

    public void setOutputVCF(boolean vcf) {
        this.VCFmode = vcf;
    }

    IndexedIdentifier columnIds = new IndexedIdentifier();
    Formatter formatter;

    /**
     * Define a new column, built with String.format and optional parameters.
     *
     * @param columnId column id format
     * @param option   options for formatting
     * @return the index of the new column
     */
    public int defineColumn(String columnId, String... option) {
        String id = String.format(columnId, option);
        final int index = columnIds.registerIdentifier(new MutableString(id));

        values = new CharSequence[columnIds.size()];
        return index;
    }

    /**
     * Define column atributes for this columnId and option. The shortId, description and columnId are formated with
     * the options to yield each column information.
     *
     * @param shortId
     * @param numValues
     * @param type
     * @param description
     * @param columnId
     * @param option
     */
    public void defineColumnAttributes(String shortId,
                                       int numValues, ColumnType type,
                                       String description,
                                       String columnId,
                                       String... option) {
        String id = String.format(columnId, option);
        shortId = String.format(shortId, option);
        description = String.format(description, option);

        final int index = columnIds.registerIdentifier(new MutableString(id));

        indexTypes.put(index, type);
        indexDescriptions.put(index, description);
        indexNumValues.put(index, numValues);
        indexShortIds.put(index, shortId);
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


    /**
     * Define new columns, built with String.format and optional parameters.
     * Columns are defined in groups of ids with different options for the same id following each other.
     *
     * @param options set of String.format ids parameters (exactly one per id)
     * @param type    Type of the new column
     * @param ids     Set of identifier format to create new column identifiers.
     */
    public void defineColumnSet(String[] options, ColumnType type, String... ids) {

        for (String columnId : ids) {   //define groups of columns, grouped by option format
            for (String option : options) {
                String id = String.format(columnId, option);
                int index = columnIds.registerIdentifier(new MutableString(id));
                indexTypes.put(index, type);
                values = new CharSequence[columnIds.size()];
            }
        }
    }

    /**
     * Define column attributes for a group of column identifiers. Description and short ids are taken to be columnId.
     * @param numValues
     * @param type
     * @param options
     * @param columnIds
     */
    public void defineColumnAttributes(int numValues, ColumnType type, String[] options, String... columnIds) {
        for (String columnId : columnIds) {   //define groups of columns, grouped by option format
            for (String option : options) {
                defineColumnAttributes(columnId, numValues, type, columnId, columnId, option);
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
        if (VCFmode) {
            outWriter.printf("##fileformat=VCFv4.1%n" +
                    "##Goby=%s%n", VersionUtils.getImplementationVersion(GobyDriver.class));
        }
        final IntCollection intCollection = columnIds.values();

        int max = -1;
        for (int columnIndex : intCollection) {
            max = Math.max(columnIndex, max);
        }
        int index = 0;
        DoubleIndexedIdentifier reverse = new DoubleIndexedIdentifier(columnIds);

        MutableString tsvHeaderLine = new MutableString();

        for (int columnIndex = 0; columnIndex <= max; columnIndex++) {
            final MutableString columnId = reverse.getId(columnIndex);
            if (VCFmode) {
                assert indexTypes.get(columnIndex) != null : "type must be defined for VCF column. Call defineColumnAttribute before writting the header.";
                outWriter.printf("##%s=<ID=%s,Number=%d,Type=%s,Description=\"%s\">%n",
                        columnId,
                        indexShortIds.get(columnIndex),
                        indexNumValues.get(columnIndex),
                        indexTypes.get(columnIndex),
                        indexDescriptions.get(columnIndex));
            }

            tsvHeaderLine.append(columnId);
            if (index != max) {
                tsvHeaderLine.append("\t");
            } else {
                tsvHeaderLine.append("\n");
            }
            ++index;
        }
        if (VCFmode) {
            outWriter.print("#");
        }
        outWriter.println(tsvHeaderLine);
        outWriter.flush();
    }

}
