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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.readers.vcf.*;
import edu.cornell.med.icb.goby.stats.*;
import edu.cornell.med.icb.io.TSVReader;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.IOUtils;

import java.io.*;
import java.util.Collections;

/**
 * Combines tab delimited datasets and performs FDR adjustment on a set of P-value columns. Lines will always be ordered
 * in the output in the same order that the lines are read from the input. However, since each line is independent, this
 * mode garantees that sorting the output by a  identifier column (unique for each line) will yield the same output
 * irrespective of the order in which the input files are presented to the mode. The FDR adjustment is done with all the
 * P-value kept in memory, but only the P-values. The data files are scanned a second time to read other columns and
 * produce the combined output.
 *
 * @author Fabien Campagne
 * @since Goby 1.9
 */
public class FalseDiscoveryRateMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "fdr";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Combines tab delimited datasets and performs FDR adjustment on a set of P-value columns.";

    /**
     * The output file.
     */
    private String outputFilename;

    /**
     * The basename of the compact alignment.
     */
    private String[] basenames;

    private double qValueThreshold;
    private String[] inputFiles;
    private String[] selectedPValueColumns;
    private ObjectArraySet<String> adjustedColumnIds = new ObjectArraySet<String>();
    private boolean vcf;


    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFiles = jsapResult.getStringArray("input");
        outputFilename = jsapResult.getString("output");
        qValueThreshold = jsapResult.getDouble("q-threshold");
        selectedPValueColumns = jsapResult.getStringArray("column");
        vcf = jsapResult.getBoolean("vcf");
        return this;
    }


    /**
     * Combine tab delimited files and adjusts P-values for multiple testing.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        Writer stream = null;
        try {
            stream = outputFilename == null ? new OutputStreamWriter(System.out)
                    : new FileWriter(outputFilename);
            DifferentialExpressionResults data = new DifferentialExpressionResults();
            ObjectList<String> columnIdList = vcf ? getVCFColumns(inputFiles) : getTSVColumns(inputFiles);
            for (String col : columnIdList) {
                System.out.println("column: " + col);
            }
            DifferentialExpressionCalculator deCalculator = new DifferentialExpressionCalculator();

            if (vcf) {
                loadVCF(inputFiles, data, deCalculator, columnIdList);
            } else {
                loadTSV(inputFiles, data, deCalculator, columnIdList);
            }
            //    data.write(new PrintWriter(System.out), '\t', deCalculator);
            BenjaminiHochbergAdjustment fdr = new BenjaminiHochbergAdjustment();
            for (String column : selectedPValueColumns) {

                fdr.adjust(data, column.toLowerCase());
            }
            for (int i = 0; i < data.getNumberOfStatistics(); i++) {
                MutableString statisticId = data.getStatisticIdForIndex(i);
                final String statId = statisticId.toString();
                if (!contains(selectedPValueColumns, statId)) {
                    adjustedColumnIds.add(statId);
                }
            }
            Collections.sort(data, new ElementIndexComparator());
            if (vcf) {
                combineVCF(inputFiles, data, deCalculator, columnIdList, stream);

            } else {
                combineTSV(inputFiles, data, deCalculator, columnIdList, stream);
            }


        } finally {
            if (outputFilename != null) {
                IOUtils.closeQuietly(stream);
            }
        }
    }


    private ObjectList<String> getVCFColumns(String[] inputFiles) {
        return new ObjectArrayList<String>();
    }


    private boolean contains(String[] selectedPValueColumns, String statId) {
        for (String s : selectedPValueColumns) {
            if (s.equalsIgnoreCase(statId)) return true;
        }
        return false;
    }

    private void loadVCF(String[] inputFiles, DifferentialExpressionResults data,
                         DifferentialExpressionCalculator deCalculator, ObjectList<String> columnIdList) throws FileNotFoundException {
        int elementIndex = 0;
        for (String filename : inputFiles) {
            VCFParser parser = new VCFParser(new FileReader(filename));
            try {
                parser.readHeader();
                IntSet selectedInfoFieldGlobalIndices = new IntArraySet();
                // find the global field indices for the INFO fields we need to load:

                for (String selectedFieldName : selectedPValueColumns) {
                    ColumnInfo infoColumn = parser.getColumns().find("INFO");
                    if (infoColumn == null) {
                        System.err.printf("Could not find INFO colum in file %s.",
                                filename);
                        System.exit(1);
                    }
                    ColumnField selectedField = infoColumn.fields.find(selectedFieldName);
                    if (selectedField == null) {
                        System.err.printf("Could not find selection colum %s in file %s.", selectedFieldName,
                                filename);
                        System.exit(1);
                    }
                    final String statName = selectedFieldName.toLowerCase();
                    if (!data.isStatisticDefined(new MutableString(statName))) {
                        data.declareStatistic(statName);
                    }
                    selectedInfoFieldGlobalIndices.add(selectedField.globalFieldIndex);
                }
                while (parser.hasNextDataLine()) {
                    parser.next();
                    // prepare elementId and stat info:

                    final String elementId = Integer.toString(elementIndex++);
                    deCalculator.defineElement(elementId);
                    int index = 0;
                    DifferentialExpressionInfo info = new DifferentialExpressionInfo(elementId);
                    info.statistics().size(selectedInfoFieldGlobalIndices.size());

                    for (int globalFieldIndex : selectedInfoFieldGlobalIndices) {

                        info.statistics().set(index++, Double.parseDouble(parser.getStringFieldValue(globalFieldIndex)));
                    }
                    data.add(info);
                }
            } catch (VCFParser.SyntaxException
                    e) {
                System.err.println("An error occured parsing VCF file " + filename);
                e.printStackTrace();
                System.exit(1);
            }
        }

    }

    private void loadTSV(String[] inputFiles, DifferentialExpressionResults data, DifferentialExpressionCalculator deCalculator,
                         ObjectList<String> columnIdList) throws IOException {
        int elementIndex = 0;
        for (String filename : inputFiles) {
            System.out.println("Loading P-values from " + filename);
            TSVReader reader = new TSVReader(new FileReader(filename), '\t');
            try {
                String firstColumn = columnIdList.get(0);
                reader.setCommentPrefix(firstColumn);

                int columnIndex = 0;
                IntSet doubleColumnIndices = new IntArraySet();

                for (String column : columnIdList) {
                    for (String selectedColumn : selectedPValueColumns) {
                        if (column.equalsIgnoreCase(selectedColumn)) {

                            final String statName = column.toLowerCase();
                            if (!data.isStatisticDefined(new MutableString(statName))) {
                                data.declareStatistic(statName);
                            }
                            doubleColumnIndices.add(columnIndex);

                        }
                    }
                    columnIndex++;
                }

                while (reader.hasNext()) {


                    if (!reader.isCommentLine()) {
                        reader.next();
                        final String elementId = Integer.toString(elementIndex++);
                        deCalculator.defineElement(elementId);
                        DifferentialExpressionInfo info = new DifferentialExpressionInfo(elementId);
                        info.statistics().size(doubleColumnIndices.size());
                        int index = 0;

                        for (int j = 0; j < reader.numTokens(); j++) {
                            if (doubleColumnIndices.contains(j)) {
                                info.statistics().set(index++, reader.getDouble());
                            } else {
                                reader.getString();
                            }

                        }
                        data.add(info);

                    } else {
                        reader.skip();
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(1);
            } finally {
                reader.close();
            }
        }

    }

    private void combineVCF(String[] inputFiles, DifferentialExpressionResults data,
                            DifferentialExpressionCalculator deCalculator,
                            ObjectList<String> columnIdList, Writer writer) throws IOException {

        Columns columns = new Columns();
        ObjectArrayList<String> sampleIdList = new ObjectArrayList();


        for (String filename : inputFiles) {
            VCFParser parser = new VCFParser(filename);
            try {
                try {
                    parser.readHeader();

                } catch (VCFParser.SyntaxException e) {
                    throw new InternalError("this syntax error should have been caught in the first pass.");
                }
                Columns fileColumns = parser.getColumns();
                for (ColumnInfo col : fileColumns) {
                    if (columns.find(col.getColumnName()) == null) {
                        columns.add(col);
                        if (col.useFormat) {
                            final String sampleId = col.getColumnName();
                            if (!sampleIdList.contains(sampleId)) {

                                sampleIdList.add(sampleId);
                            }
                        }
                    }
                }
            } finally {
                parser.close();
            }
        }
        VCFWriter vcfWriter = new VCFWriter(writer);
        vcfWriter.defineSchema(columns);
        vcfWriter.defineSamples(sampleIdList.toArray(new String[sampleIdList.size()]));

        Int2IntMap statIndexToInfoFieldIndex = new Int2IntOpenHashMap();

        // add adjusted columns:
        ColumnInfo infoColumn = columns.find("INFO");
        int statIndex = 0;
        for (String fieldName : selectedPValueColumns) {
            int newColFieldIndex = vcfWriter.defineField("INFO", fieldName + "_q", 1, ColumnType.Float,
                    String.format("Benjamini Hochberg FDR adjusted for column %s.", fieldName));
            statIndexToInfoFieldIndex.put(statIndex++, newColFieldIndex);
        }

        vcfWriter.writeHeader();
        int elementIndex = 0;
        for (String filename : inputFiles) {
            VCFParser parser = new VCFParser(filename);
            try {

                try {
                    parser.readHeader();
                } catch (VCFParser.SyntaxException e) {
                    throw new InternalError("this syntax error should have been caught in the first pass.");
                }
                final int chromosomeFieldIndex = columns.find("CHROM").getField("VALUE").globalFieldIndex;
                final int positionFieldIndex = columns.find("POS").getField("VALUE").globalFieldIndex;
                final int idFieldIndex = columns.find("ID").getField("VALUE").globalFieldIndex;
                final int refFieldIndex = columns.find("REF").getField("VALUE").globalFieldIndex;
                final int altFieldIndex = columns.find("ALT").getField("VALUE").globalFieldIndex;
                final int qualFieldIndex = columns.find("QUAL").getField("VALUE").globalFieldIndex;
                final int filterFieldIndex = columns.find("FILTER").getField("VALUE").globalFieldIndex;

                final IntSet infoFieldGlobalIndices = new IntArraySet();

                for (ColumnField infoField : parser.getColumns().find("INFO").fields) {
                    infoFieldGlobalIndices.add(infoField.globalFieldIndex);
                }

                final IntSet formatFieldGlobalIndices = new IntArraySet();
                Int2IntMap globalIndexToSampleIndex = new Int2IntOpenHashMap();
                int sampleIndex = 0;

                for (ColumnInfo col : columns) {
                    if (col.useFormat) {

                        for (ColumnField field : col.fields) {
                            globalIndexToSampleIndex.put(field.globalFieldIndex, sampleIndex++);
                            formatFieldGlobalIndices.add(field.globalFieldIndex);
                        }
                    }
                }

                int infoFieldIndex = 0;
                int formatFieldCount = 0;
                int previousSampleIndex = -1;

                while (parser.hasNextDataLine()) {

                    final String elementId = Integer.toString(elementIndex);
                    boolean keepThisLine = false;
                    for (String adjustedColumn : adjustedColumnIds) {
                        int adjustedColumnIndex = data.getStatisticIndex(adjustedColumn);
                        final DoubleArrayList list = data.get(elementIndex).statistics();

                        double adjustedPValue = list.get(adjustedColumnIndex);
                        if (adjustedPValue < qValueThreshold) {
                            keepThisLine = true;

                        }
                    }
                    if (adjustedColumnIds.size() == 0) {
                        // no selection column, keep all lines.
                        keepThisLine = true;
                    }
                    if (keepThisLine) {
                        final DifferentialExpressionInfo info = data.get(elementIndex);
                        assert info.getElementId().equals(elementId) : " elementId must match";
                        // transfer previous columsn and fields:
                        infoFieldIndex = 0;
                        sampleIndex = 0;
                        formatFieldCount = 0;
                        previousSampleIndex = -1;

                        String format = parser.getStringColumnValue(columns.find("FORMAT").columnIndex);
                        final String[] formatTokens = format.split(":");
                        int numFormatFields = columns.find("FORMAT").fields.size();
                        sampleIndex = 0;
                        int formatFieldIndex = 0;
                        for (int globalFieldIndex = 0; globalFieldIndex < parser.countAllFields(); globalFieldIndex++) {
                            String value = parser.getStringFieldValue(globalFieldIndex);

                            if (globalFieldIndex == chromosomeFieldIndex) {
                                vcfWriter.setChromosome(value);
                            } else if (globalFieldIndex == positionFieldIndex) {
                                vcfWriter.setPosition(Integer.parseInt(value));
                            } else if (globalFieldIndex == idFieldIndex) {
                                vcfWriter.setId(value);
                            } else if (globalFieldIndex == refFieldIndex) {
                                vcfWriter.setReferenceAllele(value);
                            } else if (globalFieldIndex == altFieldIndex) {
                                vcfWriter.setAlternateAllele(value);
                            } else if (globalFieldIndex == qualFieldIndex) {
                                vcfWriter.setQual(value);
                            } else if (globalFieldIndex == filterFieldIndex) {
                                vcfWriter.setFilter(value);
                            }
                            if (infoFieldGlobalIndices.contains(globalFieldIndex)) {
                                vcfWriter.setInfo(infoFieldIndex++, value);
                            }

                            if (formatFieldGlobalIndices.contains(globalFieldIndex)) {



                             /*   System.out.printf("Set sampleValue formatIndex: %d sampleIndex: %d value: %s%n", formatFieldCount, sampleIndex,
                                        value);
                                System.out.flush();
                               */ 

                                if (formatFieldIndex < formatTokens.length) {
                                    vcfWriter.setSampleValue(formatTokens[formatFieldIndex], sampleIndex, value);
                                }

                                formatFieldCount++;
                                if (value.length() != 0) {
                                    formatFieldIndex++;
                                }
                                if (formatFieldCount == numFormatFields) {
                                    formatFieldIndex = 0;
                                    formatFieldCount = 0;
                                    sampleIndex++;
                                }
                            }
                        }
                        // add new INFO field values (the adjusted p-values):
                        for (statIndex = 0; statIndex < adjustedColumnIds.size(); statIndex++) {
                            double newColValue = info.statistics().get(statIndex);
                            int infoStatFieldIndex = statIndexToInfoFieldIndex.get(statIndex);
                            vcfWriter.setInfo(infoStatFieldIndex, Double.toString(newColValue));
                        }

                        // This is a line we keep, write it:
                        vcfWriter.writeRecord();

                        elementIndex++;
                    }
                    parser.next();
                }
            } finally {
                parser.close();
            }
        }
    }

    private void combineTSV(String[] inputFiles, DifferentialExpressionResults data,
                            DifferentialExpressionCalculator deCalculator, ObjectList<String> columnIdList,
                            Writer out) throws IOException {

        PrintWriter printer = new PrintWriter(out);
        int elementIndex = 0;
        // write the TSV header first:

        for (String column : columnIdList) {
            printer.print(column);
            printer.write('\t');
        }

        boolean first = true;
        for (String column : adjustedColumnIds) {
            if (!first) {
                printer.write('\t');
            }
            printer.print(column);
            first = false;

        }
        printer.println();
        // end of TSV header generation

        for (String filename : inputFiles) {
            System.out.printf("Writing combined output (processing %s)%n", filename);
            TSVReader reader = new TSVReader(new FileReader(filename), '\t');

            try {
                String firstColumn = columnIdList.get(0);
                reader.setCommentPrefix(firstColumn);

                int columnIndex = 0;
                IntSet doubleColumnIndices = new IntArraySet();

                for (String column : columnIdList) {
                    for (String selectedColumn : selectedPValueColumns) {
                        if (column.equalsIgnoreCase(selectedColumn)) {

                            final String statName = column.toLowerCase();

                            if (!data.isStatisticDefined(new MutableString(statName))) {
                                data.declareStatistic(statName);
                            }
                            doubleColumnIndices.add(columnIndex);
                        }
                    }
                    columnIndex++;
                }

                while (reader.hasNext()) {


                    if (!reader.isCommentLine()) {
                        reader.next();
                        final String elementId = Integer.toString(elementIndex);
                        boolean keepThisLine = false;
                        for (String adjustedColumn : adjustedColumnIds) {
                            int adjustedColumnIndex = data.getStatisticIndex(adjustedColumn);
                            final DoubleArrayList list = data.get(elementIndex).statistics();

                            double adjustedPValue = list.get(adjustedColumnIndex);
                            if (adjustedPValue < qValueThreshold) {
                                keepThisLine = true;

                            }
                        }
                        if (!keepThisLine) {
                            //     System.out.println("skipping elementId since the adjusted P-values do not make the q-value threshold." + elementId);
                        }
                        if (keepThisLine) {
                            int index = 0;
                            final DifferentialExpressionInfo info = data.get(elementIndex);
                            assert info.getElementId().equals(elementId) : " elementId must match";
                            for (int j = 0; j < reader.numTokens(); j++) {
                                if (doubleColumnIndices.contains(j)) {
                                    reader.getString();

                                    printer.print(info.statistics().get(index));
                                    printer.print('\t');
                                    index++;
                                } else {
                                    printer.print(reader.getString());
                                    printer.print('\t');
                                }

                            }
                            first = true;
                            for (String adjustedColumn : adjustedColumnIds) {
                                int adjustedColumnIndex = data.getStatisticIndex(adjustedColumn);
                                final DoubleArrayList list = data.get(elementIndex).statistics();

                                if (!first) {
                                    printer.write('\t');
                                }
                                printer.print(list.get(adjustedColumnIndex));
                                first = false;
                            }
                            printer.printf("%n");
                        }

                        elementIndex++;

                    } else {
                        reader.skip();
                    }
                }
            } catch (IOException
                    e) {
                e.printStackTrace();
                System.exit(1);
            } finally {
                reader.close();
            }
        }

        printer.flush();
    }

    private ObjectList<String> getTSVColumns
            (String[] inputFiles) throws IOException {
        ObjectArrayList<String> columns = new ObjectArrayList<String>();
        for (String filename : inputFiles) {
            TSVReader reader = new TSVReader(new FileReader(filename), '\t');
            try {
                if (reader.hasNext()) {
                    reader.next();

                    // int numColumns = reader.numTokens();
                    //     for (int i = 0; i < numColumns; i++) {
                    for (int i = 0; i < reader.numTokens(); i++) {
                        final String s = reader.getString().trim();
                        String cols[] = s.split("[\t]");
                        for (String c : cols) {
                            // System.out.println("col: "+c);
                            c = c.trim();
                            if (!columns.contains(c)) columns.add(c);
                        }
                    }
                }
            } finally {
                reader.close();
            }
        }
        return columns;
    }

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */

    public static void main
            (
                    final String[] args) throws JSAPException, IOException {
        new FalseDiscoveryRateMode().configure(args).execute();
    }
}