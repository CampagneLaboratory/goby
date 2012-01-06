/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.readers.vcf;

import edu.cornell.med.icb.goby.modes.TabToColumnInfoMode;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;
import net.sf.samtools.util.BlockCompressedInputStream;
import org.apache.commons.io.IOUtils;

import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;

/**
 * Parser for files in the <a href="http://vcftools.sourceforge.net/specs.html">Variant Call Format</a>, or in plain TSV format.
 * This parser will read either a VCF file with meta-data information defined about columns, or will read plain TSV files.
 * Two version of TSV files are supported. That which starts with a header line without # or with # as first character.
 * This reader supports a Group attribute in field declarations. When the VCF 4 supports the following declaration line,
 * which declares ID, Number, Type and Description attributes:
 * <pre>##INFO=&lt;ID=AF1,Number=1,Type=Float,Description="Max-likelihood ..."&gt;</pre>
 * This parser can additionally read a Group attribute, such as in:
 * <pre>##INFO=&lt;ID=AF1,Number=1,Type=Float,Group=LIKELIHOODS,Description="Max-likelihood ..."&gt;</pre>
 *
 * @author Fabien Campagne
 *         Date: Mar 26, 2011
 *         Time: 3:01:47 PM
 */
public class VCFParser implements Closeable {
    private Reader input;
    private Columns columns = new Columns();
    private boolean hasNextDataLine;
    private int numberOfColumns;
    private int[] columnStarts;
    private int[] columnEnds;
    private MutableString line;
    private char columnSeparatorCharacter = '\t';
    private int[] fieldStarts;
    private int[] fieldEnds;
    private char fieldSeparatorCharacter = ';';
    private char formatFieldSeparatorCharacter = ':';
    private int numberOfFields;
    private int globalFieldIndex;
    private int formatColumnIndex;
    private int lineLength;
    private int globalColumnIndex;
    /**
     * Variable TSV is true if we determined the file is tab delimited.
     */
    private boolean TSV = true;
    /**
     * An array with dimensions numAllFields that stores the permutation from the global field index
     * to the observed field index (taking into account absence or presence of Flag attributes, and the fact that
     * fields that occur in any order on each line in a column).
     */
    private int[] fieldPermutation;
    /**
     * Sorts columns in increasing columnIndex order.
     */
    public static final Comparator<ColumnInfo> COLUMN_SORT = new Comparator<ColumnInfo>() {
        @Override
        public int compare(final ColumnInfo c1, final ColumnInfo c2) {
            return c1.columnIndex - c2.columnIndex;
        }
    };
    /**
     * Sorts fields in increasing globalFieldIndex order.
     */
    public static final Comparator<ColumnField> FIELD_SORT = new Comparator<ColumnField>() {
        @Override
        public int compare(final ColumnField c1, final ColumnField c2) {
            return c1.globalFieldIndex - c2.globalFieldIndex;
        }
    };
    private final ObjectArrayList<ColumnInfo> columnList = new ObjectArrayList<ColumnInfo>();
    private final ObjectArrayList<ColumnField> fieldList = new ObjectArrayList<ColumnField>();
    private ColumnInfo formatColumn;
    private boolean headerLineNotParsed = true;
    private boolean headerParsed;
    private File inputFile;

    /**
     * Map of column name to column type for TSV files. When reading a TSV file this is populated by
     * readTsvColumnTypes() during parseHeaderLine() and used to set the correct types based on the
     * data found in the TSV file.
     */
    private Map<String, ColumnType> tsvColumnNameToTypeMap;

    /**
     * If we should cache the TSV column types to a .colinfo file when parsing TSV files. This must be
     * set BEFORE calling readHeader() if you want it to be used.
     */
    private boolean cacheTsvColumnTypes = true;

    /**
     * When scanning a TSV file to determine column types, this is the number of lines that will be checked.
     * Set to <= 0 to scan the entire file. This must be set before calling readHeader() for the value to be used.
     */
    private int tsvLinesToScanForColumnType = 10000;
    /**
     * Indicate that the field permutation has already been computed. When this field is true, the permutation
     * of the first line computed is reused for subsequent lines.
     */
    private boolean computedFieldPermutation;
    /**
     * When this field is true, the permutation of the first line computed is reused for subsequent lines.
     */
    private boolean cacheFieldPermutation;

    /**
     * Constructs a VCF parser.
     *
     * @param file Input to parse
     */
    public VCFParser(final Reader file) {
        this.input = file;
    }

    /**
     * Constructs a VCF parser.  When the filename ends in .gz, this method attempts to decompress
     * the file on the fly, using BlockCompressedInputStream (from samtools). BlockCompressedInputStream
     * is used preferentially to GZipInputStream to avoid truncating bgzip input files produced with bgzip.
     * http://biostar.stackexchange.com/questions/6112/how-to-decompress-1000genomes-bgzip-compressed-files-using-java
     *
     * @param filename Input to parse
     * @throws java.io.IOException when an error occurs.
     */
    public VCFParser(final String filename) throws IOException {
        inputFile = new File(filename);
        input = filename.endsWith(".gz") ?
                new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(filename))) :
                new FileReader(filename);
    }

    /**
     * If we should cache the TSV column types to a .colinfo file when parsing TSV files. This must be set
     * before calling readHeader() for the value to be used.
     *
     * @param cacheTsvColumnTypes new cacheTsvColumnTypes value
     */
    public void setCacheTsvColumnTypes(final boolean cacheTsvColumnTypes) {
        this.cacheTsvColumnTypes = cacheTsvColumnTypes;
    }

    /**
     * If we should cache the TSV column types to a .colinfo file when parsing TSV files.
     *
     * @return cacheTsvColumnTypes value
     */
    public boolean isCacheTsvColumnTypes() {
        return cacheTsvColumnTypes;
    }

    /**
     * When scanning a TSV file to determine column types, this is the number of lines that will be checked.
     * Set to <= 0 to scan the entire file. This must be set before calling readHeader() for the value to be used.
     *
     * @param tsvLinesToScanForColumnType the new value of tsvLinesToScanForColumnType
     */
    public void setTsvLinesToScanForColumnType(final int tsvLinesToScanForColumnType) {
        this.tsvLinesToScanForColumnType = tsvLinesToScanForColumnType;
    }

    /**
     * When scanning a TSV file to determine column types, this is the number of lines that will be checked.
     *
     * @return the value of tsvLinesToScanForColumnType
     */
    public int getTsvLinesToScanForColumnType() {
        return tsvLinesToScanForColumnType;
    }

    /**
     * If parsing a TSV file, this can be called to retrieve the column types map. The column types map
     * .colinfo cache file will be created the first time this is called and read on subsequent calls.
     *
     * @throws IOException if error reading / creating .colinfo cache file
     */
    public void readTsvColumnTypes() throws IOException {
        // TODO: Currently this will not work right if the TSV file header line DOES start with "#"
        // TODO: which could be a problem, although current TSV files produced by DiffExp's do not do this.
        if (tsvColumnNameToTypeMap == null && inputFile != null && inputFile.exists()) {
            final TabToColumnInfoMode reader = new TabToColumnInfoMode();
            reader.addInputFile(inputFile);
            reader.setCreateCache(cacheTsvColumnTypes);
            reader.setNumberOfLinesToProcess(tsvLinesToScanForColumnType);
            reader.setReadFromCache(true);
            reader.execute();
            tsvColumnNameToTypeMap = reader.getDetailsAtIndex(0);
        }
    }

    /**
     * Return the number of columns in the file. This method can be called after the header has been read to obtain
     * the number of columns in the file.
     *
     * @return The number of columns in the file
     */
    public int getNumberOfColumns() {
        return numberOfColumns;
    }

    /**
     * Return the columns in the file. This method can be called after the header has been read.
     *
     * @return The columns declared in the file
     */
    public Columns getColumns() {
        return columns;
    }

    public ColumnType getColumnType(final int columnIndex) {
        for (final ColumnInfo col : columns) {
            if (col.columnIndex == columnIndex) {
                if (col.fields.size() == 1) {
                    return ((ColumnField) col.fields.toArray()[0]).type;

                } else {
                    break;
                }

            }
        }

        return ColumnType.String;

    }

    /**
     * Return the type of the specified field.
     *
     * @param globalFieldIndex global field index, from zero to countAllFields()-1
     * @return type of the specified field.
     */
    public ColumnType getFieldType(final int globalFieldIndex) {
        return fieldList.get(globalFieldIndex).type;
    }

    /**
     * Return the number of values in the specified field.
     *
     * @param globalFieldIndex global field index, from zero to countAllFields()-1
     * @return the number of values contained in this field.
     */
    public int getFieldNumValues(final int globalFieldIndex) {
        return fieldList.get(globalFieldIndex).numberOfValues;
    }

    public String getColumnName(final int columnIndex) {
        for (final ColumnInfo col : columns) {
            if (col.columnIndex == columnIndex) {
                return col.columnName;

            }
        }

        return null;

    }

    private LineIterator lineIterator;

    public boolean hasNextDataLine() {
        if (hasNextDataLine) {
            return true;
        }

        hasNextDataLine = lineIterator.hasNext();
        if (hasNextDataLine) {
            line = lineIterator.next();
            if (!TSV) {
                parseCurrentLine();

            } else {
                parseTSVLine();
            }

        }
        return hasNextDataLine;
    }


    public void next() {
        if (!hasNextDataLine) {


            throw new IllegalArgumentException("Next can be called only after hasNext has returned true.");
        }
        hasNextDataLine = false;
    }

    /**
     * Returns a column value as a CharSequence. Faster than returning a String.
     *
     * @param columnIndex index of the field on a line of input.
     * @return a column value as a CharSequence.
     */
    public CharSequence getColumnValue(final int columnIndex) {
        if (hasNextDataLine) {

            return line.subSequence(columnStarts[columnIndex], columnEnds[columnIndex]);

        } else return null;
    }

    /**
     * Returns the total number of fields across all columns.
     *
     * @return the sum of the number of fields in each column.
     */
    public int countAllFields() {
        int n = 0;
        for (final ColumnInfo column : columns) {

            n += column.fields.size();
        }

        return n;
    }

    /**
     * Returns the value of the field.
     * The field is identified by a global index that runs from zero (inclusive) to countAllFields() (exclusive).
     *
     * @param globalFieldIndex a global index that runs from zero to countAllFields()
     * @return Value of this field.
     */
    public CharSequence getFieldValue(final int globalFieldIndex) {
        if (hasNextDataLine) {

            final int lineFieldIndex = fieldPermutation[globalFieldIndex];
            if (lineFieldIndex == -1) {
                // missing field in this row;
                return "";
            }
            final int start = fieldStarts[lineFieldIndex];
            final int end = fieldEnds[lineFieldIndex];

            assert (start >= 0 && end <= lineLength) :
                    String.format("position indices must be within line boundaries start: %d end: %d length: %d", start, end, lineLength);
            return line.subSequence(start, end);


        } else return null;
    }

    /**
     * Returns the value of a field.
     * The field is identified by a global index that runs from zero (inclusive) to countAllFields() (exclusive).
     *
     * @param globalFieldIndex a global index that runs from zero to countAllFields()
     * @return Value of this field.
     */
    public String getStringFieldValue(final int globalFieldIndex) {

        final CharSequence value = getFieldValue(globalFieldIndex);
        return value == null ? null : value.toString();
    }

    /**
     * Returns a column value as a String.
     *
     * @param columnIndex index of the field on a line of input.
     * @return a column value as a String.
     */
    public String getStringColumnValue(final int columnIndex) {
        return getColumnValue(columnIndex).toString();
    }

    private Int2ObjectMap<String> fieldIndexToName;

    /**
     * Return the field name, in the format:
     * <LI>For columns with multiple fields: &lt;column-name&gt;[&lt;field-id&gt;]
     * <LI>For columns with a single field: &lt;column-name&gt;[&lt;field-id&gt;]
     *
     * @param globalFieldIndex index of the field across columns.
     * @return the field name
     */
    public String getFieldName(final int globalFieldIndex) {
        return fieldIndexToName.get(globalFieldIndex);
    }

    private FastBufferedReader bufferedReader;

    /**
     * Read the header of this file. Headers in the VCF format are supported, as well as TSV single header lines (with or
     * without first character #.
     *
     * @throws SyntaxException When the syntax of the VCF file is incorrect.
     */
    public void readHeader() throws SyntaxException {
        if (headerParsed) {
            return;
        }
        headerParsed = true;
        globalFieldIndex = 0;
        fieldIndexToName = new Int2ObjectOpenHashMap<String>();
        bufferedReader = new FastBufferedReader(input);
        lineIterator = new LineIterator(bufferedReader);
        int lineNumber = 1;
        try {
            while (lineIterator.hasNext()) {
                line = lineIterator.next();
                if (line.startsWith("##")) {
                    TSV = false;
                }
                if (!line.startsWith("#")) {
                    if (TSV && lineNumber == 1 && headerLineNotParsed) {

                        // assume the file is TSV and starts directly with the header line. Parse lineIterator here.
                        parseHeaderLine(new MutableString("#" + line));

                    } else {
                        // We are seeing an actual line of data. Prepare for parsing:
                        if (!TSV) {
                            parseCurrentLine();
                        } else {
                            parseTSVLine();
                        }
                        hasNextDataLine = true;
                    }
                    break;
                }
                if (line.startsWith("##")) {
                    TSV = false;
                    processMetaInfoLine(line);
                } else if (line.startsWith("#")) {
                    parseHeaderLine(line);
                }
                lineNumber++;
            }
        } catch (net.sf.samtools.FileTruncatedException e) {
            line = null;
            hasNextDataLine = false;
            // install a dummy line iteratory that always returns false when hasNext is called.
            lineIterator = new LineIterator(bufferedReader) {
                @Override
                public boolean hasNext() {
                    return false;
                }

                @Override
                public MutableString next() {
                    return null;
                }
            };

        }
    }

    private void parseTSVLine() {
        Arrays.fill(columnStarts, 0);
        Arrays.fill(columnEnds, 0);
        Arrays.fill(fieldStarts, 0);
        Arrays.fill(fieldEnds, 0);
        int columnIndex = 0;
        lineLength = line.length();
        for (int i = 0; i < lineLength; i++) {

            final char c = line.charAt(i);
            if (c == '\t') {
                final String columnName = columnList.get(columnIndex).columnName;

                columnEnds[columnIndex] = i;
                if (columnIndex + 1 < columnStarts.length) {

                    columnStarts[columnIndex + 1] = i + 1;
                }
                fieldPermutation[columnIndex] = columnIndex;
                ++columnIndex;
            }
        }

        fieldPermutation[columnEnds.length - 1] = columnEnds.length - 1;
        columnEnds[columnEnds.length - 1] = lineLength;
        columnStarts[columnEnds.length - 1] = columnEnds[columnEnds.length - 2] + 1;
        System.arraycopy(columnEnds, 0, fieldEnds, 0, columnEnds.length);
        System.arraycopy(columnStarts, 0, fieldStarts, 0, columnStarts.length);

    }

    final IntArrayList previousColumnFieldIndices = new IntArrayList();

    private void parseCurrentLine() {
        Arrays.fill(columnStarts, 0);
        Arrays.fill(columnEnds, 0);
        Arrays.fill(fieldStarts, 0);
        Arrays.fill(fieldEnds, 0);
        columnStarts[0] = 0;
        int columnIndex = 0;
        int fieldIndex = 0;
        lineLength = line.length();
        final int[] lineFieldIndexToColumnIndex = new int[numberOfFields];
        Arrays.fill(lineFieldIndexToColumnIndex, -1);
        previousColumnFieldIndices.clear();
        // determine the position of column and field delimiters:
        final char[] chrs = line.toCharArray();
        for (int i = 0; i < lineLength; i++) {
            final char c = chrs[i];
            if (c == columnSeparatorCharacter) {

                fieldPermutation[columnIndex] = columnIndex;
                columnEnds[columnIndex] = i;

                if (columnIndex + 1 < numberOfColumns) {
                    columnStarts[columnIndex + 1] = i + 1;
                }
                //lineFieldIndexToColumnIndex[columnIndex] = columnIndex;
            }
            if (c == columnSeparatorCharacter ||
                    c == fieldSeparatorCharacter ||
                    (columnIndex >= formatColumnIndex &&
                            c == formatFieldSeparatorCharacter)) {

                if (TSV) {

                    // there are no fields, only columns, the field separators do not apply

                    fieldEnds[columnIndex] = columnEnds[columnIndex];
                    fieldStarts[columnIndex] = columnStarts[columnIndex];
                    fieldIndex = columnIndex;
                    lineFieldIndexToColumnIndex[fieldIndex] = columnIndex;
                } else {

                    fieldEnds[fieldIndex] = i;

                    if (fieldIndex + 1 < numberOfFields) {
                        fieldStarts[fieldIndex + 1] = i + 1;
                    }

                    previousColumnFieldIndices.add(fieldIndex);
                    fieldIndex++;
                }
            }
            if (c == columnSeparatorCharacter) {
                if (TSV) {
                    lineFieldIndexToColumnIndex[fieldIndex] = columnIndex;
                }
                push(columnIndex, lineFieldIndexToColumnIndex, previousColumnFieldIndices);
                columnIndex++;
            }

        }
        int numberOfFieldsOnLine = Math.min(fieldIndex, fieldEnds.length-1);
        int numberOfColumnsOnLine = Math.min(columnIndex,columnEnds.length-1);
        columnStarts[0] = 0;
        columnEnds[numberOfColumnsOnLine - (TSV ? 1 : 0)] = line.length();
        fieldStarts[0] = 0;
        fieldEnds[numberOfFieldsOnLine - (TSV ? 1 : 0)] = line.length();
        previousColumnFieldIndices.add(fieldIndex);
        push(columnIndex, lineFieldIndexToColumnIndex, previousColumnFieldIndices);

        if (cacheFieldPermutation && computedFieldPermutation) return;
        Arrays.fill(fieldPermutation, -1);
        for (ColumnInfo c : columns) {
            c.formatIndex = 0;
        }


        // determine the fieldPermutation for each possible field:
        for (int lineFieldIndex = 0; lineFieldIndex <= numberOfFieldsOnLine; lineFieldIndex++) {

            final int start = fieldStarts[lineFieldIndex];
            final int end = fieldEnds[lineFieldIndex];

            final int cIndex = lineFieldIndexToColumnIndex[lineFieldIndex];
            if (cIndex >= columnList.size()) {
                break;
            }
            final ColumnInfo column = columnList.get(cIndex);

            int colMinGlobalFieldIndex = Integer.MAX_VALUE;
            int colMaxGlobalFieldIndex = Integer.MIN_VALUE;
            final ColumnFields fields = column.fields;
            fields.rebuildList();

            for (int fi = 0; fi < fields.size(); ++fi) {

                final ColumnField f = fields.get(fi);
                colMinGlobalFieldIndex = Math.min(colMinGlobalFieldIndex, f.globalFieldIndex);
                colMaxGlobalFieldIndex = Math.max(colMaxGlobalFieldIndex, f.globalFieldIndex);

            }

            final int formatColumnIndex = TSV ? -1 : formatColumn.columnIndex;
            final int startFormatColumn = TSV ? 0 : columnStarts[formatColumnIndex];
            final int endFormatColumn = TSV ? 0 : columnEnds[formatColumnIndex];


            final String[] formatTokens = split(line, formatFieldSeparatorCharacter, startFormatColumn, endFormatColumn);

            for (int fi = 0; fi < fields.size(); ++fi) {

                final ColumnField f = fields.get(fi);

                if (fieldPermutation[f.globalFieldIndex] != -1) {
                    // already assigned.
                    continue;
                }
                if (colMaxGlobalFieldIndex == colMinGlobalFieldIndex) {
                    // This column has only one field.
                    fieldPermutation[f.globalFieldIndex] = lineFieldIndex;
                    break;
                } else {
                    // find the column field f whose id matches the character span we are looking at :
                    int j = start;
                    final String id = f.id;
                    int matchLength = 0;
                    for (int i = 0; i < id.length(); i++) {
                        if (j >= end) {
                            // reached end of field, not this field.
                            break;
                        }

                        final char linechar = line.charAt(j);

                        if (id.charAt(i) != linechar) {
                            // found mimatch with field id, not this field.
                            matchLength = -1;
                            break;
                        }
                        matchLength++;
                        j++;
                    }

                    if (matchLength == id.length() && line.charAt(j) == '=' ||
                            (j == end && f.type == ColumnType.Flag)) {
                        // found the correct field.
                        /*  System.out.printf("Assigning global %s %d -> %d for field %s%n",
                                f.id, globalFieldIndex, lineFieldIndex, line.subSequence(start, end));
                        */
                        fieldPermutation[f.globalFieldIndex] = lineFieldIndex;
                        if (f.type != ColumnType.Flag) {
                            fieldStarts[lineFieldIndex] += f.id.length() + 1; // remove id= from value;
                            //fieldStarts[lineFieldIndex]=Math.min(fieldStarts[lineFieldIndex],fieldEnds[lineFieldIndex]);
                        }
                        break;
                    } else {

                        if (column.useFormat && column.formatIndex < formatTokens.length) {

                            if (f.id.equals(formatTokens[column.formatIndex])) {
                                /*    System.out.printf("Assigning FORMAT global %s %d -> %d for field %s%n",
                                          f.id, f.globalFieldIndex, lineFieldIndex, line.subSequence(start, end));
                                */
                                fieldPermutation[f.globalFieldIndex] = lineFieldIndex;
                                column.formatIndex++;
                                break;
                            }
                        }
                    }
                }
            }
        }
        computedFieldPermutation = true;
    }

    String[] formatSplit = null;

    private String[] split(final MutableString line, final char formatFieldSeparatorCharacter,
                           final int startFormatColumn, final int endFormatColumn) {
        if (formatSplit != null) {
            return formatSplit;
        } else {
            final MutableString formatSpan = line.substring(startFormatColumn, endFormatColumn);

            int fieldCount = 0;

            formatSpan.append(formatFieldSeparatorCharacter);
            final int length = formatSpan.length();
            for (int i = 0; i < length; i++) {
                if (formatSpan.charAt(i) == formatFieldSeparatorCharacter) {
                    ++fieldCount;
                }
            }
            final String[] result = new String[fieldCount];
            final MutableString value = new MutableString();
            int last = 0;
            int j = 0;
            for (int i = 0; i < length; i++) {
                if (formatSpan.charAt(i) == formatFieldSeparatorCharacter && i > last) {
                    value.append(formatSpan.substring(last, i));
                    last = i + 1;
                    result[j] = value.toString();
                    value.setLength(0);
                    ++j;
                }
            }
            formatSplit = result;
            return result;
        }
    }

    //     System.out.println("ned");


    private void push(final int columnIndex, final int[] lineFieldIndexToColumnIndex, final IntArrayList previousColumnFieldIndices) {
        //    System.out.println("---");
        final int size = previousColumnFieldIndices.size();
        for (int i = 0; i < size; ++i) {
            final int fIndex = previousColumnFieldIndices.getInt(i);
            /*        System.out.printf("field %s gfi:%d belongs to column %d %s%n ",
           line.substring(fieldStarts[fIndex], fieldEnds[fIndex]),
           fIndex,
           columnIndex, columnList.get(columnIndex).columnName);*/
            lineFieldIndexToColumnIndex[fIndex] = columnIndex;
        }
        previousColumnFieldIndices.clear();
    }


    private void parseHeaderLine(MutableString line) {
        if (TSV) {
            // Attempt to determine the column types using TabToColumnInfoMode.
            try {
                readTsvColumnTypes();
            } catch (IOException e) {
                System.err.println("Could not determine column info from tsv file " + e.getMessage());
            }
        }

        headerLineNotParsed = false;
        // System.out.printf("header line:%s%n", line);
        // drop the #
        line = line.substring(1);
        final String[] columnNames = line.toString().split("[\\t]");

        for (final String columnName : columnNames) {

            defineFixedColumn(columnName);
            if (!columns.hasColumnName(columnName)) {
                final ColumnInfo formatColumn = columns.find("FORMAT");

                // copy the fields of the FORMAT column for each sample:
                final ColumnField[] fields;
                if (formatColumn != null) {
                    fields = new ColumnField[formatColumn.fields.size()];
                    int i = 0;
                    for (final ColumnField field : formatColumn.fields) {
                        fields[i] = new ColumnField(field.id, field.numberOfValues, field.type, field.description);
                        fields[i].globalFieldIndex = -1;
                        i++;
                    }
                } else {
                    final ColumnType columnType;
                    if (tsvColumnNameToTypeMap != null) {
                        columnType = tsvColumnNameToTypeMap.get(columnName) == null ?
                                ColumnType.String :
                                tsvColumnNameToTypeMap.get(columnName);
                    } else {
                        columnType = ColumnType.String;
                    }
                    fields = new ColumnField[]{new ColumnField("VALUE", 1, columnType, "")};
                    fields[0].globalFieldIndex = -1;
                }
                final ColumnInfo newCol = new ColumnInfo(columnName, fields);
                newCol.useFormat = true;
                columns.add(newCol);
            }

        }
        formatColumn = columns.find("FORMAT");
        // columns.remove(formatColumn);
        // columnList.remove(formatColumn);

        for (final ColumnInfo column : columns) {
            if (column.columnIndex == -1) {
                column.columnIndex = globalColumnIndex++;
            }
            for (final ColumnField field : column.fields) {

                if (field.globalFieldIndex == -1) {
                    field.globalFieldIndex = globalFieldIndex++;
                }
                final String name;
                if (column.fields.size() == 1) {
                    name = column.columnName;

                } else {
                    name = String.format("%s[%s]", column.columnName, field.id);
                }
                fieldIndexToName.put(field.globalFieldIndex, name);
            }
        }
        formatColumnIndex = TSV ? -1 : formatColumn.columnIndex;

        numberOfColumns = globalColumnIndex;
        columnStarts = new int[numberOfColumns];
        columnEnds = new int[numberOfColumns];
        numberOfFields = globalFieldIndex;
        fieldStarts = new int[numberOfFields];
        fieldEnds = new int[numberOfFields];
        fieldPermutation = new int[numberOfFields];

        columnList.addAll(columns);
        Collections.sort(columnList, COLUMN_SORT);
        for (final ColumnInfo column : columnList) {
            fieldList.addAll(column.fields);
        }
    }

    private void defineFixedColumn(final String columnName) {
        for (final ColumnInfo fixed : fixedColumns) {
            if (fixed.columnName.equals(columnName) && !columns.hasColumnName(columnName)) {
                fixed.columnIndex = globalColumnIndex++;
                for (final ColumnField field : fixed.fields) {
                    field.globalFieldIndex = globalFieldIndex++;
                }
                columns.add(fixed);
                return;
            }
        }


    }

    public static ColumnInfo[] fixedColumn() {

        // returns a deep copy of fixed columns so that this class is not affected by possible changes to the returned
        // value.
        final ColumnInfo[] copy = new ColumnInfo[fixedColumns.length];
        int i = 0;

        for (final ColumnInfo fixedColumn : fixedColumns) {
            copy[i++] = fixedColumn.copy();
        }

        return copy;

    }

    final static private ColumnInfo fixedColumns[] = new ColumnInfo[]{
            new ColumnInfo("CHROM", new ColumnField("VALUE", 1, ColumnType.String, "The reference position, with the 1st base having position 1. " +
                    "Positions are sorted numerically, in increasing order, within each reference sequence CHROM.")),
            new ColumnInfo("POS", new ColumnField("VALUE", 1, ColumnType.Integer, "he reference position, with the 1st base having position 1. " +
                    "Positions are sorted numerically, in increasing order, within each reference sequence CHROM.")),
            new ColumnInfo("ID", new ColumnField("VALUE", 1, ColumnType.String, "ID semi-colon separated list of unique identifiers where available. " +
                    "If this is a dbSNP variant it is encouraged to use the rs number(s). " +
                    "No identifier should be present in more than one data record. " +
                    "If there is no identifier available, then the missing value should be used.")),
            new ColumnInfo("REF", new ColumnField("VALUE", 1, ColumnType.String, "Reference base(s): Each base must be one of A,C,G,T,N. " +
                    "Bases should be in uppercase. Multiple bases are permitted. " +
                    "The value in the POS field refers to the position of the first base in the String. " +
                    "For InDels, the reference String must include the base before the event " +
                    "(which must be reflected in the POS field).")),

            new ColumnInfo("ALT", new ColumnField("VALUE", 1, ColumnType.String,
                    "Comma separated list of alternate non-reference alleles called on at least one of the samples. " +
                            "Options are base Strings made up of the bases A,C,G,T,N, or " +
                            "an angle-bracketed ID String (\"<ID>\"). " +
                            "If there are no alternative alleles, then the missing value should be used. " +
                            "Bases should be in uppercase. (Alphanumeric String; no whitespace, commas, " +
                            "or angle-brackets are permitted in the ID String itself).")),
            new ColumnInfo("QUAL", new ColumnField("VALUE", 1, ColumnType.Float,
                    "Phred-scaled quality score for the assertion made in ALT. i.e. give -10log_10 prob(call in ALT is wrong). " +
                            "If ALT is \".\" (no variant) then this is -10log_10 p(variant), " +
                            "and if ALT is not \".\" this is -10log_10 p(no variant). " +
                            "High QUAL scores indicate high confidence calls. " +
                            "Although traditionally people use integer phred scores, this field is permitted to be " +
                            "a floating point to enable higher resolution for low confidence calls if desired.")),
            new ColumnInfo("FILTER", new ColumnField("VALUE", 1, ColumnType.String,
                    "Filter: PASS if this position has passed all filters, i.e. a call is made at this position. " +
                            "Otherwise, if the site has not passed all filters, a semicolon-separated list of codes " +
                            "for filters that fail. e.g. \"10;s50\" might indicate that at this site the quality is " +
                            "below 10 and the number of samples with data is below 50%% of the total number of samples. " +
                            "\"0\" is reserved and should not be used as a filter String. " +
                            "If filters have not been applied, then this field should be set to the missing value.")),
            new ColumnInfo("INFO", new ColumnField("VALUE", 1, ColumnType.String,
                    "Additional information: INFO fields are encoded as a semicolon-separated series of short keys " +
                            "with optional values in the format: <key>=<data>[,data]. Arbitrary keys are permitted, " +
                            "although some sub-fields are reserved.")),
    };


    private void processMetaInfoLine(final MutableString line) throws SyntaxException {
        final int start = 2;
        final int end = line.indexOf('=');
        final String columnName = line.substring(start, end).toString();

        final MutableString restOfLine = line.substring(end + 1);
        processMetaInfo(columnName, restOfLine);
    }


    private void processMetaInfo(final String columnName, final MutableString infoDefinition) throws SyntaxException {
        if ("fileformat".equals(columnName) ||
                "samtoolsVersion".equals(columnName)) {
            return;
        }
        if (!infoDefinition.startsWith("<") && !infoDefinition.endsWith(">")) {
            // is this a syntax error? test-data/vcf/tricky.vcf would trigger errors
            return;
        }
        final ColumnInfo info;
        if (columns.hasColumnName(columnName)) {
            info = columns.find(columnName);
        } else {

            info = new ColumnInfo();
            info.columnName = columnName;
            columns.add(info);
        }

        final ColumnField field = new ColumnField();
        final MutableString insideBrackets = infoDefinition.substring(1, infoDefinition.length() - 1);

        try {
            final String[] tokens = insideBrackets.toString().split("(,N)|(,T)|(,D)|(,G)");
            for (final String token : tokens) {
                final String[] kv = new String[2];
                final int firstEqualIndex = token.indexOf('=');
                kv[0] = token.substring(0, firstEqualIndex);
                kv[1] = token.substring(firstEqualIndex + 1);

                if ("ID".equals(kv[0])) {
                    field.id = kv[1];
                } else if ("umber".equals(kv[0])) {

                    String num = kv[1];
                    if (".".equals(num)) {
                        // . indicates any number of values. Use 2
                        num = "-1";
                    }
                    field.numberOfValues = Integer.parseInt(num);
                } else if ("ype".equals(kv[0])) {
                    field.type = ColumnType.valueOf(kv[1]);
                } else if ("roup".equals(kv[0])) {
                    field.group = kv[1];
                } else if ("escription".equals(kv[0])) {
                    if (kv[1].startsWith("\"") && kv[1].endsWith("\"")) {
                        kv[1] = kv[1].substring(1, kv[1].length() - 1);
                    }
                    field.description = kv[1];
                } else {
                    throw new SyntaxException(infoDefinition);
                }

            }
        } catch (NumberFormatException e) {
            throw new SyntaxException(infoDefinition);
        }
        // System.out.println("adding " + field);
        // do not set the global field index on a meta-info field yet. We will do this after fixed columns have been added.
        info.addField(field);
    }

    /**
     * Return a global field index, or -1 if the column or field id are not declared.
     *
     * @param columnName name of column.
     * @param fieldId    Identifier for field in column.
     * @return a global field index, or -1 if the column or field id are not declared.
     */
    public int getGlobalFieldIndex(final String columnName, final String fieldId) {
        final ColumnInfo column = columns.find(columnName);
        if (column == null) {
            return -1;
        }
        final ColumnField columnField = column.fields.find(fieldId);
        if (columnField == null) {
            return -1;
        }
        return columnField.globalFieldIndex;
    }

    /**
     * Releases the IO resources held by this parser.
     *
     * @throws IOException
     */
    @Override
    public void close() throws IOException {
        if (bufferedReader != null) {
            IOUtils.closeQuietly(bufferedReader);

        }
        IOUtils.closeQuietly(input);

    }

    /**
     * Return the sample names, or more specifically, the names of column that use the FORMAT column.
     *
     * @return
     */
    public String[] getColumnNamesUsingFormat() {
        final ObjectArrayList<String> columnNamesUsingFormat = new ObjectArrayList<String>();
        for (final ColumnInfo info : columns) {
            if (info.useFormat) {
                columnNamesUsingFormat.add(info.columnName);
            }
        }
        return columnNamesUsingFormat.toArray(new String[columnNamesUsingFormat.size()]);
    }

    public void setCacheFieldPermutation(boolean cacheFieldPermutation) {
        this.cacheFieldPermutation = cacheFieldPermutation;
    }


    public class SyntaxException extends Exception {
        public SyntaxException(final MutableString line) {
            super(line.toString());
        }
    }


}
