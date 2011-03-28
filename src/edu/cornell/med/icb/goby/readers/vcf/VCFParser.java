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

import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;

import java.io.Reader;

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
public class VCFParser {
    private Reader input;
    private Columns columns = new Columns();
    private boolean hasNextDataLine;
    private int numberOfColumns;
    private int[] columnStarts;
    private int[] columnEnds;
    private MutableString line;
    private char separatorCharacter = '\t';

    /**
     * Constructs a VCF parser.
     *
     * @param file Input to parse
     */
    public VCFParser(Reader file) {
        this.input = file;
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

    private void parseCurrentLine() {
        columnStarts[0] = 0;
        int fieldIndex = 0;
        for (int i = 0; i < line.length(); i++) {
            if (line.charAt(i) == separatorCharacter) {


                columnEnds[fieldIndex] = i;

                if (fieldIndex + 1 < numberOfColumns) {
                    columnStarts[fieldIndex + 1] = i + 1;
                }
                fieldIndex++;
            }
        }
    }

    public void next() {
        if (!hasNextDataLine)
            throw new IllegalArgumentException("Next can be called only after hasNext has returned true.");
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
     * Returns a column value as a String.
     *
     * @param columnIndex index of the field on a line of input.
     * @return a column value as a String.
     */
    public String getStringColumnValue(int columnIndex) {
        return getColumnValue(columnIndex).toString();
    }

    /**
     * Read the header of this file. Headers in the VCF format are supported, as well as TSV single header lines (with or
     * without first character #.
     * @throws SyntaxException When the syntax of the VCF file is incorrect.
     */
    public void readHeader() throws SyntaxException {
        lineIterator = new LineIterator(new FastBufferedReader(input));
        int lineNumber = 1;
        while (lineIterator.hasNext()) {
            line = lineIterator.next();
            if (!line.startsWith("#")) {
                if (lineNumber == 1) {
                    // assume the file is TSV and starts directly with the header line. Parse lineIterator here.
                    processHeaderLine(new MutableString("#" + line));
                } else {
                    // We are seeing an actual line of data. Prepare for parsing:
                    parseCurrentLine();
                    hasNextDataLine = true;
                }
                break;
            }
            if (line.startsWith("##")) {
                processMetaInfoLine(line);
            } else if (line.startsWith("#")) {
                processHeaderLine(line);
            }
            lineNumber++;
        }
    }

    private void processHeaderLine(MutableString line) {
        // System.out.printf("header line:%s%n", line);
        // drop the #
        line = line.substring(1);
        String[] columnNames = line.toString().split("[\\s]");
        int columnIndex = 0;
        for (String columnName : columnNames) {
            ColumnInfo column = columns.find(columnName);
            if (column != null) {
                column.columnIndex = columnIndex;
            }
            defineFixedColumn(columnName, columnIndex);
            if (!columns.hasColumnName(columnName)) {

                ColumnInfo newCol = new ColumnInfo(columnName, new ColumnField("VALUE", -1,
                        ColumnField.ColumnType.String, "Format described for each row in the associated FORMAT column."));

                newCol.columnIndex = columnIndex;
                columns.add(newCol);
            }
            columnIndex++;
        }
        numberOfColumns = columnIndex;
        columnStarts = new int[numberOfColumns];
        columnEnds = new int[numberOfColumns];
    }

    private void defineFixedColumn(String columnName, int columnIndex) {
        for (ColumnInfo fixed : fixedColumns) {
            if (fixed.columnName.equals(columnName) && !columns.hasColumnName(columnName)) {
                fixed.columnIndex = columnIndex;
                columns.add(fixed);
                return;
            }
        }


    }


    final private ColumnInfo fixedColumns[] = new ColumnInfo[]{
            new ColumnInfo("CHROM", new ColumnField("VALUE", 1, ColumnField.ColumnType.String, "The reference position, with the 1st base having position 1. " +
                    "Positions are sorted numerically, in increasing order, within each reference sequence CHROM.")),
            new ColumnInfo("POS", new ColumnField("VALUE", 1, ColumnField.ColumnType.Integer, "he reference position, with the 1st base having position 1. " +
                    "Positions are sorted numerically, in increasing order, within each reference sequence CHROM.")),
            new ColumnInfo("ID", new ColumnField("VALUE", 1, ColumnField.ColumnType.String, "ID semi-colon separated list of unique identifiers where available. " +
                    "If this is a dbSNP variant it is encouraged to use the rs number(s). " +
                    "No identifier should be present in more than one data record. " +
                    "If there is no identifier available, then the missing value should be used.")),
            new ColumnInfo("REF", new ColumnField("VALUE", 1, ColumnField.ColumnType.String, "Reference base(s): Each base must be one of A,C,G,T,N. " +
                    "Bases should be in uppercase. Multiple bases are permitted. " +
                    "The value in the POS field refers to the position of the first base in the String. " +
                    "For InDels, the reference String must include the base before the event " +
                    "(which must be reflected in the POS field).")),

            new ColumnInfo("ALT", new ColumnField("VALUE", 1, ColumnField.ColumnType.String,
                    "Comma separated list of alternate non-reference alleles called on at least one of the samples. " +
                            "Options are base Strings made up of the bases A,C,G,T,N, or " +
                            "an angle-bracketed ID String (Ó<ID>Ó). " +
                            "If there are no alternative alleles, then the missing value should be used. " +
                            "Bases should be in uppercase. (Alphanumeric String; no whitespace, commas, " +
                            "or angle-brackets are permitted in the ID String itself).")),
            new ColumnInfo("QUAL", new ColumnField("VALUE", 1, ColumnField.ColumnType.Float,
                    "Phred-scaled quality score for the assertion made in ALT. i.e. give -10log_10 prob(call in ALT is wrong). " +
                            "If ALT is Ó.Ó (no variant) then this is -10log_10 p(variant), " +
                            "and if ALT is not Ó.Ó this is -10log_10 p(no variant). " +
                            "High QUAL scores indicate high confidence calls. " +
                            "Although traditionally people use integer phred scores, this field is permitted to be " +
                            "a floating point to enable higher resolution for low confidence calls if desired.")),
            new ColumnInfo("FILTER", new ColumnField("VALUE", 1, ColumnField.ColumnType.String,
                    "Filter: PASS if this position has passed all filters, i.e. a call is made at this position. " +
                            "Otherwise, if the site has not passed all filters, a semicolon-separated list of codes " +
                            "for filters that fail. e.g. Òq10;s50Ó might indicate that at this site the quality is " +
                            "below 10 and the number of samples with data is below 50% of the total number of samples. " +
                            "Ò0Ó is reserved and should not be used as a filter String. " +
                            "If filters have not been applied, then this field should be set to the missing value.")),
            new ColumnInfo("INFO", new ColumnField("VALUE", 1, ColumnField.ColumnType.String,
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
            throw new SyntaxException(infoDefinition);
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


        String tokens[] = insideBrackets.toString().split("(,N)|(,T)|(,D)|(,G)");
        for (String token : tokens) {
            String[] kv = new String[2];
            final int firstEqualIndex = token.indexOf("=");
            kv[0] = token.substring(0, firstEqualIndex);
            kv[1] = token.substring(firstEqualIndex + 1);

            if ("ID".equals(kv[0])) {
                field.id = kv[1];
            } else if ("umber".equals(kv[0])) {
                field.numberOfValues = Integer.parseInt(kv[1]);
            } else if ("ype".equals(kv[0])) {
                field.type = ColumnField.ColumnType.valueOf(kv[1]);
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
        // System.out.println("adding " + field);
        info.addField(field);
    }


    public class SyntaxException extends Exception {
        public SyntaxException(MutableString line) {
            super(line.toString());
        }
    }

    
}
