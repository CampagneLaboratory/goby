/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

import edu.cornell.med.icb.goby.modes.GobyDriver;
import edu.cornell.med.icb.goby.readers.vcf.*;
import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;

/**
 * Helper class to write VCF statistic result files.
 *
 * @author Fabien Campagne
 *         Date: Sep 23, 2010
 *         Time: 1:26:56 PM
 */
public class VCFWriter {
    PrintWriter outWriter;
    private Int2ObjectMap<ColumnType> indexTypes;
    private boolean VCFmode;
    private Int2ObjectMap<String> indexDescriptions;
    private Int2IntMap indexNumValues;
    private Int2ObjectMap<String> indexShortIds;

    public int chromosomeColumnIndex;
    public int positionColumnIndex;
    private CharSequence chrom;
    private int position;
    private String id;
    private int infoFieldIndex = 0;
    private CharSequence[] formatFieldIds;
    private CharSequence ref;
    private CharSequence alt;
    private CharSequence filter;

    private boolean[] formatFieldActive;
    private CharSequence[][] formatValues;
    private CharSequence[] infoIds;


    public void setId(String id) {
        this.id = id;
    }

    public void setPosition(int position) {

        this.position = position;
    }

    public void setChromosome(CharSequence chromosome) {
        this.chrom = chromosome;
    }

    private Columns columns = new Columns();
    private ObjectArrayList<ColumnInfo> columnList = new ObjectArrayList<ColumnInfo>();

    public VCFWriter(PrintWriter outWriter) {
        this.outWriter = outWriter;
        this.sampleIds=new String[0];
        columns.addAll(Arrays.asList(VCFParser.fixedColumn()));
        for (ColumnInfo c : columns) {
            for (ColumnField field : c.fields) {
                if (field.id.equals("VALUE")) {
                    c.fields.remove(field);
                }
            }
        }
        columns.add(new ColumnInfo("FORMAT"));
    }

    public void writeHeader() {

        outWriter.printf("##fileformat=VCFv4.1%n" +
                "##Goby=%s%n", VersionUtils.getImplementationVersion(GobyDriver.class));
        columnList.addAll(columns);
        Collections.sort(columnList, VCFParser.COLUMN_SORT);


        int index = 0;

        MutableString tsvHeaderLine = new MutableString();

        for (ColumnInfo c : columnList) {

            for (ColumnField field : c.fields) {
                outWriter.printf("##%s=<ID=%s,Number=%d,Type=%s,Description=\"%s\">%n",
                        c.getColumnName(),
                        field.id,
                        field.numberOfValues,
                        field.type,
                        field.description);
            }
            tsvHeaderLine.append(c.getColumnName());

            tsvHeaderLine.append("\t");

            ++index;
        }


        for (String sampleId : sampleIds) {

            tsvHeaderLine.append(sampleId);
            tsvHeaderLine.append("\t");
            ++index;

        }
        // trim last \t
        tsvHeaderLine.setLength(tsvHeaderLine.length() - 1);
        tsvHeaderLine.append("\n");
        outWriter.print("#");
        outWriter.print(tsvHeaderLine);
        outWriter.flush();

        infoValues = new CharSequence[columns.find("INFO").fields.size()];
        final int numFormatTypes = columns.find("FORMAT").fields.size();
        formatFieldIds = new CharSequence[numFormatTypes];
        formatFieldActive = new boolean[numFormatTypes];
        index = 0;
        for (ColumnField formatField : columns.find("FORMAT").fields) {
            formatFieldIds[index++] = formatField.id;
        }
        formatValues = new CharSequence[formatFieldActive.length][sampleIds.length];
        ref = "";
        alt = "";
        filter = "";
        id = "";
        chrom = "";
        final ColumnInfo info = columns.find("INFO");
        infoIds = new CharSequence[info.fields.size()];

        for (int infoFieldIndex = 0; infoFieldIndex < infoIds.length; infoFieldIndex++) {

            infoIds[infoFieldIndex] = info.
                    fields.find(infoFieldIndex).id;
        }
    }

    public void writeRecord() {
        outWriter.append(chrom);
        outWriter.append('\t');
        outWriter.append(Integer.toString(position));
        outWriter.append('\t');
        outWriter.append(id);
        outWriter.append('\t');
        outWriter.append(ref);
        outWriter.append('\t');
        outWriter.append(alt);
        outWriter.append('\t');
        outWriter.append(filter);
        outWriter.append('\t');
        int max;
        int index = 0;
        max = infoValues.length;
        for (CharSequence infoValue : infoValues) {
            outWriter.append(infoIds[index]);
            outWriter.append('=');
            outWriter.append(infoValue);

            if (++index != max) outWriter.append(';');

        }
        outWriter.append('\t');
        int formatIndex = 0;
        max = formatFieldIds.length;
        index = 0;
        for (CharSequence formatType : formatFieldIds) {

            if (formatFieldActive[formatIndex]) {
                outWriter.append(formatFieldIds[formatIndex]);
                if (++index != max) {
                    outWriter.append(':');
                }
            }
            formatIndex++;
        }
        outWriter.append('\t');
        int sampleIndex = 0;

        max = formatFieldIds.length;
        for (CharSequence value : sampleIds) {
            formatIndex = 0;
            index = 0;

            for (CharSequence formatType : formatFieldIds) {
                if (formatFieldActive[formatIndex]) {
                    outWriter.append(formatValues[formatIndex][sampleIndex]);
                    if (++index != max) outWriter.append(':');

                }
                formatIndex++;

            }
            sampleIndex++;
            outWriter.append('\t');
        }
        outWriter.println();
        Arrays.fill(formatFieldActive, false);
        for (int i = 0; i < formatFieldActive.length; i++) Arrays.fill(formatValues[i], "");
        Arrays.fill(infoValues, "");
        ref = "";
        alt = "";
        filter = "";
        id = "";
        chrom = "";

    }


    public static int COLUMN_NOT_DEFINED = -1;

    public void close() {
        outWriter.close();
    }

    CharSequence[] infoValues;


    ColumnInfo infoColumn = new ColumnInfo("INFO");


    public int defineField(String columnName, String fieldName, int numValues, ColumnType type, String description) {
        ColumnInfo c = columns.find(columnName);
        if (c == null) {
            throw new IllegalArgumentException("Could not find column " + columnName);
        }
        int maxIndex = -1;
        for (ColumnField f : c.fields) {
            maxIndex = Math.max(f.globalFieldIndex, maxIndex);
        }
        final ColumnField columnField = new ColumnField(fieldName, numValues, type, description);
        columnField.globalFieldIndex = maxIndex + 1;
        c.addField(columnField);

        return columnField.globalFieldIndex;
    }

    public void setInfo(int infoFieldIndex, CharSequence value) {
        infoValues[infoFieldIndex] = value;
    }

    public void setFormat(int formatFieldIndex, CharSequence fieldId) {
        formatFieldIds[formatFieldIndex] = fieldId;
    }

    String[] sampleIds;

    public void defineSamples(String[] samples) {
        sampleIds = samples;
        formatValues = new CharSequence[getNumFormatFields()][samples.length];
    }


    public void setSampleValue(int fieldIndex, int sampleIndex, CharSequence value) {
        formatFieldActive[fieldIndex] = true;
        formatValues[fieldIndex][sampleIndex] = value;
    }

    public int getNumFormatFields() {
        return columns.find("FORMAT").fields.size();
    }
}