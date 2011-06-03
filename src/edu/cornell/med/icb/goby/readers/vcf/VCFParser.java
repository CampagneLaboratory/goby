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

package edu.cornell.med.icb.goby.readers.vcf;

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.IOUtils;

import java.io.*;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.zip.GZIPInputStream;

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
        public int compare(ColumnInfo c1, ColumnInfo c2) {
            return c1.columnIndex - c2.columnIndex;
        }
    };
    /**
     * Sorts fields in increasing globalFieldIndex order.
     */
    public static final Comparator<ColumnField> FIELD_SORT = new Comparator<ColumnField>() {
        public int compare(ColumnField c1, ColumnField c2) {
            return c1.globalFieldIndex - c2.globalFieldIndex;
        }
    };
    private ObjectArrayList<ColumnInfo> columnList = new ObjectArrayList<ColumnInfo>();
    private ObjectArrayList<ColumnField> fieldList = new ObjectArrayList<ColumnField>();
    private ColumnInfo formatColumn;
    private boolean headerLineNotParsed = true;
    private boolean headerParsed;


    /**
     * Constructs a VCF parser.
     *
     * @param file Input to parse
     */
    public VCFParser(Reader file) {
        this.input = file;
    }

    /**
     * Constructs a VCF parser.
     *
     * @param filename Input to parse
     */
    public VCFParser(String filename) throws IOException {
        this.input = filename.endsWith(".gz") ?
                new InputStreamReader(new GZIPInputStream(new FileInputStream(filename), 100000)) :
                new FileReader(filename);
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

    public ColumnType getColumnType(int columnIndex) {
        for (ColumnInfo col : columns) {
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
    public ColumnType getFieldType(int globalFieldIndex) {
        return fieldList.get(globalFieldIndex).type;
    }

    /**
     * Return the number of values in the specified field.
     *
     * @param globalFieldIndex global field index, from zero to countAllFields()-1
     * @return the number of values contained in this field.
     */
    public int getFieldNumValues(int globalFieldIndex) {
        return fieldList.get(globalFieldIndex).numberOfValues;
    }

    public String getColumnName(int columnIndex) {
        for (ColumnInfo col : columns) {
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
            parseCurrentLine();

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
    public CharSequence getColumnValue(int columnIndex) {
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
        for (ColumnInfo column : columns) {

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
    public String getStringFieldValue(int globalFieldIndex) {

        final CharSequence value = getFieldValue(globalFieldIndex);
        return value == null ? null : value.toString();
    }

    /**
     * Returns a column value as a String.
     *
     * @param columnIndex index of the field on a line of input.
     * @return a column value as a String.
     */
    public String getStringColumnValue(int columnIndex) {
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
    public String getFieldName(int globalFieldIndex) {
        return fieldIndexToName.get(globalFieldIndex);
    }

    private FastBufferedReader bufferedReader = null;

    /**
     * Read the header of this file. Headers in the VCF format are supported, as well as TSV single header lines (with or
     * without first character #.
     *
     * @throws SyntaxException When the syntax of the VCF file is incorrect.
     */
    public void readHeader() throws SyntaxException {
        if (headerParsed) return;
        headerParsed = true;
        globalFieldIndex = 0;
        fieldIndexToName = new Int2ObjectOpenHashMap<String>();
        bufferedReader = new FastBufferedReader(input);
        lineIterator = new LineIterator(bufferedReader);
        int lineNumber = 1;
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
                    parseCurrentLine();
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
        int lineFieldIndexToColumnIndex[] = new int[numberOfFields];
        Arrays.fill(lineFieldIndexToColumnIndex, -1);
        previousColumnFieldIndices.clear();
        // determine the position of column and field delimiters:

        for (int i = 0; i < lineLength; i++) {
            final char c = line.charAt(i);
            if (c == columnSeparatorCharacter) {

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
                push(columnIndex, lineFieldIndexToColumnIndex, previousColumnFieldIndices);
                columnIndex++;
            }

        }
        int numberOfFieldsOnLine = fieldIndex;
        int numberOfColumnsOnLine = columnIndex;
        columnStarts[0] = 0;
        columnEnds[numberOfColumnsOnLine] = line.length();
        fieldStarts[0] = 0;
        fieldEnds[numberOfFieldsOnLine] = line.length();
        previousColumnFieldIndices.add(fieldIndex);
        push(columnIndex, lineFieldIndexToColumnIndex, previousColumnFieldIndices);


        Arrays.fill(fieldPermutation, -1);
        for (ColumnInfo c : columns) {
            c.formatIndex = 0;
        }


        // determine the fieldPermutation for each possible field:
        for (int lineFieldIndex = 0; lineFieldIndex <= numberOfFieldsOnLine; lineFieldIndex++) {

            int start = fieldStarts[lineFieldIndex];
            int end = fieldEnds[lineFieldIndex];

            final int cIndex = lineFieldIndexToColumnIndex[lineFieldIndex];

            ColumnInfo column = columnList.get(cIndex);

            int colMinGlobalFieldIndex = Integer.MAX_VALUE;
            int colMaxGlobalFieldIndex = Integer.MIN_VALUE;
            final ColumnFields fields = column.fields;
            fields.rebuildList();

            for (int fi = 0; fi < fields.size(); ++fi) {

                ColumnField f = fields.get(fi);
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
            String result[] = new String[fieldCount];
            MutableString value = new MutableString();
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
        //  System.out.println("---");
        for (final int fIndex : previousColumnFieldIndices.toIntArray()) {
            /*       System.out.printf("field %s gfi:%d belongs to column %d %s%n ",
           line.substring(fieldStarts[fIndex], fieldEnds[fIndex]),
           fIndex,
           columnIndex, columnList.get(columnIndex).columnName);*/
            lineFieldIndexToColumnIndex[fIndex] = columnIndex;
        }
        previousColumnFieldIndices.clear();
    }


    private void parseHeaderLine(MutableString line) {
        headerLineNotParsed = false;
        // System.out.printf("header line:%s%n", line);
        // drop the #
        line = line.substring(1);
        String[] columnNames = line.toString().split("[\\t]");

        for (String columnName : columnNames) {

            defineFixedColumn(columnName);
            if (!columns.hasColumnName(columnName)) {
                ColumnInfo formatColumn = columns.find("FORMAT");

                // copy the fields of the FORMAT column for each sample:
                ColumnField[] fields;
                if (formatColumn != null) {
                    fields = new ColumnField[formatColumn.fields.size()];
                    int i = 0;
                    for (ColumnField f : formatColumn.fields) {
                        fields[i] = (new ColumnField(f.id, f.numberOfValues,
                                f.type, f.description));
                        fields[i].globalFieldIndex = -1;
                        i++;
                    }
                } else {
                    fields = new ColumnField[]{new ColumnField("VALUE", 1, ColumnType.String, "")};
                    fields[0].globalFieldIndex = -1;
                }
                ColumnInfo newCol = new ColumnInfo(columnName, fields);

                newCol.useFormat = true;
                columns.add(newCol);
            }

        }
        formatColumn = columns.find("FORMAT");
        // columns.remove(formatColumn);
        // columnList.remove(formatColumn);

        for (ColumnInfo column : columns) {
            if (column.columnIndex == -1) {
                column.columnIndex = globalColumnIndex++;
            }
            for (ColumnField field : column.fields) {

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
        numberOfFields = TSV ? numberOfColumns : globalFieldIndex;
        fieldStarts = new int[numberOfFields];
        fieldEnds = new int[numberOfFields];
        fieldPermutation = new int[numberOfFields];

        columnList.addAll(columns);
        Collections.sort(columnList, COLUMN_SORT);
        for (ColumnInfo column : columnList) {
            fieldList.addAll(column.fields);
        }
    }

    private void defineFixedColumn(String columnName) {
        for (ColumnInfo fixed : fixedColumns) {
            if (fixed.columnName.equals(columnName) && !columns.hasColumnName(columnName)) {
                fixed.columnIndex = globalColumnIndex++;
                for (ColumnField field : fixed.fields) {
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
        ColumnInfo copy[] = new ColumnInfo[fixedColumns.length];
        int i = 0;

        for (ColumnInfo c : fixedColumns) {
            copy[i++] = c.copy();
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


    private void processMetaInfoLine(MutableString line) throws SyntaxException {
        int start = 2;
        int end = line.indexOf('=');
        String columnName = line.substring(start, end).toString();

        final MutableString restOfLine = line.substring(end + 1);
        processMetaInfo(columnName, restOfLine);
    }


    private void processMetaInfo(String columnName, MutableString infoDefinition) throws SyntaxException {
        if ("fileformat".equals(columnName) ||
                "samtoolsVersion".equals(columnName)) {
            return;
        }
        if (!infoDefinition.startsWith("<") && !infoDefinition.endsWith(">")) {
            // is this a syntax error? test-data/vcf/tricky.vcf would trigger errors
            return;
        }
        ColumnInfo info;
        if (columns.hasColumnName(columnName)) {
            info = columns.find(columnName);
        } else {

            info = new ColumnInfo();
            info.columnName = columnName.toString();
            columns.add(info);
        }

        ColumnField field = new ColumnField();
        final MutableString insideBrackets = infoDefinition.substring(1, infoDefinition.length() - 1);

        try {
            String tokens[] = insideBrackets.toString().split("(,N)|(,T)|(,D)|(,G)");
            for (String token : tokens) {
                String[] kv = new String[2];
                final int firstEqualIndex = token.indexOf("=");
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
    public int getGlobalFieldIndex(String columnName, String fieldId) {
        final ColumnInfo column = columns.find(columnName);
        if (column == null) return -1;
        final ColumnField columnField = column.fields.find(fieldId);
        if (columnField == null) return -1;
        return columnField.globalFieldIndex;
    }

    /**
     * Releases the IO resources held by this parser.
     *
     * @throws IOException
     */
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
        ObjectArrayList<String> columnNamesUsingFormat = new ObjectArrayList<String>();
        for (ColumnInfo info : columns) {
            if (info.useFormat) {
                columnNamesUsingFormat.add(info.columnName);
            }
        }
        return columnNamesUsingFormat.toArray(new String[columnNamesUsingFormat.size()]);
    }


    public class SyntaxException extends Exception {
        public SyntaxException(MutableString line) {
            super(line.toString());
        }
    }


}
