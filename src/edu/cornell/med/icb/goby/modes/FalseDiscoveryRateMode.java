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
import edu.cornell.med.icb.goby.stats.*;
import edu.cornell.med.icb.io.TSVReader;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
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
 *
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

        return this;
    }


    /**
     * Suggests slices to process a large alignment file in parallel.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        PrintStream stream = null;
        try {
            stream = outputFilename == null ? System.out
                    : new PrintStream(new FileOutputStream(outputFilename));
            DifferentialExpressionResults data = new DifferentialExpressionResults();
            ObjectList<String> columnIdList = getColumns(inputFiles);
            for (String col : columnIdList) {
                System.out.println("column: " + col);
  }
            DifferentialExpressionCalculator deCalculator = new DifferentialExpressionCalculator();

            load(inputFiles, data, deCalculator, columnIdList);

            data.write(new PrintWriter(System.out), '\t', deCalculator);
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
            combine(inputFiles, data, deCalculator, columnIdList, stream);


        } finally {
            if (stream != System.out) {
                IOUtils.closeQuietly(stream);
            }
        }
    }

    private boolean contains(String[] selectedPValueColumns, String statId) {
        for (String s : selectedPValueColumns) {
            if (s.equalsIgnoreCase(statId)) return true;
        }
        return false;
    }

    private void load(String[] inputFiles, DifferentialExpressionResults data, DifferentialExpressionCalculator deCalculator, ObjectList<String> columnIdList) throws IOException {
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

    private void combine(String[] inputFiles, DifferentialExpressionResults data,
                         DifferentialExpressionCalculator deCalculator, ObjectList<String> columnIdList,
                         PrintStream out) throws IOException {
        int elementIndex = 0;
        for (String column : columnIdList) {
            out.print(column);
            out.write('\t');
        }
        boolean first = true;
        for (String column : adjustedColumnIds) {
            if (!first) {
                out.write('\t');
            }
            out.print(column);
            first = false;

        }
        out.println();
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

                        int index = 0;
                        final DifferentialExpressionInfo info = data.get(elementIndex);
                        assert info.getElementId().equals(elementId) : " elementId must match";
                        for (int j = 0; j < reader.numTokens(); j++) {
                            if (doubleColumnIndices.contains(j)) {
                                reader.getString();

                                out.print(info.statistics().get(index));
                                out.print('\t');
                                index++;
                            } else {
                                out.print(reader.getString());
                                out.print('\t');
                            }

                        }
                        first = true;
                        for (String adjustedColumn : adjustedColumnIds) {
                            int adjustedColumnIndex = data.getStatisticIndex(adjustedColumn);
                            final DoubleArrayList list = data.get(elementIndex).statistics();

                            if (!first) {
                                out.write('\t');
                            }
                            out.print(list.get(adjustedColumnIndex));
                            first = false;
                        }
                        elementIndex++;
                        out.printf("%n");
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
        out.flush();
    }

    private ObjectList<String> getColumns
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