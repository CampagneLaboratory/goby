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
import java.io.Writer;
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
    private int position = -1;
    private String id;
    private CharSequence[] formatFieldIds;
    private MutableString ref;
    private MutableString alt;
    private CharSequence filter;

    private boolean[] formatFieldActive;
    private CharSequence[][] formatValues;
    private CharSequence[] infoIds;
    ObjectArrayList<String> refAlleles;
    private ObjectArrayList<String> altAlleles;
    private char genotypeDelimiterCharacter;

    public VCFWriter(Writer writer) {
        this(new PrintWriter(writer));
    }

    /**
     * Indicate whether the genotypes should be recorded as phased (true) or unphased (false).
     * Default value at construction of the writer is unphased.
     *
     * @param state True when the genotypes are phased.
     */
    public void setGenotypesPhased(boolean state) {
        genotypeDelimiterCharacter = state ? '|' : '/';
    }

    public void setId(String id) {
        this.id = id;
    }

    public void setPosition(int position) {

        this.position = position;
    }

    public void setChromosome(CharSequence chromosome) {
        this.chrom = chromosome;
    }

    public void setReferenceAllele(String allele) {
        refAlleles.clear();
        refAlleles.add(allele);
    }

    public void addAlternateAllele(String allele) {
        if (!altAlleles.contains(allele)) {
            altAlleles.add(allele);
        }
    }

    private Columns columns = new Columns();
    private ObjectArrayList<ColumnInfo> columnList = new ObjectArrayList<ColumnInfo>();

    /**
     * Contruct a VCFWriter.
     *
     * @param outWriter Where the output will be written.
     */
    public VCFWriter(PrintWriter outWriter) {
        this.outWriter = outWriter;
        this.sampleIds = new String[0];
        columns.addAll(Arrays.asList(VCFParser.fixedColumn()));
        for (ColumnInfo c : columns) {
            for (ColumnField field : c.fields) {
                if (field.id.equals("VALUE")) {
                    c.fields.remove(field);
                }
            }
        }
        columns.add(new ColumnInfo("FORMAT"));
        refAlleles = new ObjectArrayList<String>();
        altAlleles = new ObjectArrayList<String>();
        ref = new MutableString();
        alt = new MutableString();
        setGenotypesPhased(false);
    }

    /**
     * Write the VCF header.
     */
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
        ref.setLength(0);
        alt.setLength(0);
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

    /**
     * Append a record. Call the various setter before writing a record.
     */
    public void writeRecord() {
        outWriter.append(chrom);
        outWriter.append('\t');
        outWriter.append(position == -1 ? "" : Integer.toString(position));
        outWriter.append('\t');
        outWriter.append(id);
        outWriter.append('\t');
        outWriter.append(constructAlleleString(refAlleles));
        outWriter.append('\t');
        outWriter.append(constructAlleleString(altAlleles));
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

        filter = "";
        id = "";
        chrom = "";
        altAlleles.clear();
        refAlleles.clear();
        position = -1;
    }

    MutableString buffer = new MutableString();

    private MutableString constructAlleleString(ObjectArrayList<String> refAlleles) {
        buffer.setLength(0);
        int max = refAlleles.size();
        int index = 0;
        for (String allele : refAlleles) {
            buffer.append(allele);
            if (++index != max) buffer.append(',');
        }
        return buffer;
    }

    private MutableString codedGenotypeBuffer = new MutableString();

    /**
     * Encode the list of alleles as a VCF genotype.
     *
     * @param alleles String of the form A/B/C where A,B C are allele that are part of a genotype.
     * @return coded VCF genotype.
     * @see #codeGenotype(String[] alleles)
     */
    public MutableString codeGenotype(String alleles) {
        return codeGenotype(alleles.split("/"));
    }

    /**
     * Encode the list of alleles as a VCF genotype. VCF genotypes are coded against the REF and ALT alleles.
     * A VCF genotype has the form 0/2/3 where each int encodes the index of the allele that participates to the
     * genotype. REF ALT aleleles are considered in the order they appear and given increasing indices (starting at
     * zero with the REF allele). The coded genotype 0/0 represents a genotype with two reference alleles.
     *
     * @param alleles List of alleles included in the genotype.
     * @return coded VCF genotype.
     */
    public MutableString codeGenotype(String[] alleles) {
        codedGenotypeBuffer.setLength(0);
        boolean alleleFound = false;
        for (String allele : alleles) {
            int alleleIndex = 0;
            for (String ref : refAlleles) {
                if (ref.equals(allele)) {
                    codedGenotypeBuffer.append(Integer.toString(alleleIndex));
                    codedGenotypeBuffer.append(genotypeDelimiterCharacter);
                    alleleFound = true;
                    continue;
                }
                alleleIndex++;
            }
            for (String alt : altAlleles) {
                if (alt.equals(allele)) {
                    codedGenotypeBuffer.append(Integer.toString(alleleIndex));
                    codedGenotypeBuffer.append(genotypeDelimiterCharacter);
                    alleleFound = true;
                }
                alleleIndex++;
            }
            if (!alleleFound) {
                throw new IllegalArgumentException(String.format("Allele %s was not found in REF or ALT", allele));
            }
        }

        final int length = codedGenotypeBuffer.length();
        if (length > 0) {
            codedGenotypeBuffer.setLength(length - 1);
        }
        return codedGenotypeBuffer.copy();
    }

    public static int COLUMN_NOT_DEFINED = -1;

    /**
     * Close the writer.
     */
    public void close() {
        outWriter.close();
    }

    CharSequence[] infoValues;


    ColumnInfo infoColumn = new ColumnInfo("INFO");

    /**
     * Define a VCF field for a column.
     *
     * @param columnName  Name of an existing column.
     * @param fieldName   Name of the new field.
     * @param numValues   Number of values in this field.
     * @param type        Type of data for individual values of the field.
     * @param description Description of data represented by the field.
     * @return index of the field in the column.
     */
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

    /**
     * Set the value of an INFO column field for the current record. The value will be written when writeRecord is executed.
     *
     * @param infoFieldIndex Index returned by definedField("INFO",...)
     * @param value          Value of the field.
     */
    public void setInfo(int infoFieldIndex, CharSequence value) {
        infoValues[infoFieldIndex] = value;
    }


    private String[] sampleIds;

    /**
     * Get the list of sample column identifiers.
     *
     * @return list of sample column identifiers.
     */
    public String[] getSamples() {
        return sampleIds;
    }

    /**
     * Define sample identifiers. VCF stores information for each sample according to information stored in the
     * FORMAT column in each line.
     * se
     *
     * @param samples Identifiers of the samples
     */
    public void defineSamples(String[] samples) {
        sampleIds = samples;
        formatValues = new CharSequence[getNumFormatFields()][samples.length];
    }

    /**
     * Set a value of a sample column. The sampleIndex identifies the sample in the getSampleIds()  array.
     *
     * @param formatFieldIndex  Index of a FORMAT field created with defineField("FORMAT,...)
     * @param sampleIndex Index of the sample
     * @param value       Value to set the field to for the current record.
     */
    public void setSampleValue(int formatFieldIndex, int sampleIndex, CharSequence value) {
        formatFieldActive[formatFieldIndex] = true;
        formatValues[formatFieldIndex][sampleIndex] = value;
    }

    /**
     * Get the total number of possible format fields in the FORMAT column.
     *
     * @return number of FORMAT fields.
     */
    public int getNumFormatFields() {
        return columns.find("FORMAT").fields.size();
    }

    /**
     * Define the VCF schema from Columns assembled externally.
     *
     * @param columns
     */
    public void defineSchema(Columns columns) {
        columnList.clear();
        this.columns = columns;
    }

    public void setAlternateAllele(String value) {
        String[] alleles = value.split(",");
        for (String allele : alleles) {
            addAlternateAllele(allele);
        }
    }

    public void setFilter(String value) {
        this.filter=value;
    }
}