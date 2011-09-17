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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.readers.vcf.*;
import edu.cornell.med.icb.goby.stats.InformativeColumns;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.zip.GZIPOutputStream;

/**
 * Converts a VCF file to tab delimited format.
 *
 * @author Fabien Campagne
 * @since Goby 1.9.8
 */
public class VcfToTabMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "vcf-to-tab";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Converts a VCF file to tab delimited format.";

    /**
     * The output file.
     */
    private String outputFilename;


    private String[] inputFiles;
    private String[] selectedColumns;
    private final ObjectArraySet<String> selectedColumnIds = new ObjectArraySet<String>();

    private static final Logger LOG = Logger.getLogger(VcfToTabMode.class);

    public VcfToTabMode() {
    }


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
        selectedColumns = jsapResult.getStringArray("column");

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
                    : outputFilename.endsWith(".gz") ?
                    new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFilename))) :
                    new FileWriter(outputFilename);
            // start with an array of size 1M. This improves loading time for large datasets.

            for (String filename : inputFiles) {
                System.out.printf("Converting %s%n", filename);
                VCFParser parser = new VCFParser(filename);
                try {
                    parser.readHeader();
                    IntSet selectedInfoFieldGlobalIndices = new IntArraySet();
                    ColumnInfo infoColumn = parser.getColumns().find("INFO");
                    ColumnInfo formatColumn = parser.getColumns().find("FORMAT");
                    if (selectedColumns.length == 0) {

                        // add all INFO fields:

                        for (ColumnField field : infoColumn.fields) {
                            selectedColumnIds.add(field.id);
                        }
                        /*  // add all FORMAT fields:
for (ColumnField field : formatColumn.fields) {
selectedColumnIds.add(field.id);
}                                                                                 */
                        selectedColumns = selectedColumnIds.toArray(new String[selectedColumnIds.size()]);
                    }
                    // find the global field indices for the INFO fields we need to load:

                    for (String selectedFieldName : selectedColumns) {

                        if (infoColumn == null) {
                            System.err.printf("Could not find INFO column in file %s.",
                                    filename);
                            System.exit(1);
                        }
                        ColumnField selectedField = findColumnField(selectedFieldName, infoColumn, formatColumn);
                        if (selectedField == null) {
                            System.err.printf("Could not find selection column %s in file %s.", selectedFieldName,
                                    filename);
                            System.exit(1);
                        }

                        stream.write(selectedFieldName);
                        stream.write("\t");
                        selectedInfoFieldGlobalIndices.add(selectedField.globalFieldIndex);
                    }
                    stream.write("\n");
                    ProgressLogger pg = new ProgressLogger(LOG);

                    pg.priority = org.apache.log4j.Level.INFO;
                    pg.itemsName = "line";
                    pg.displayFreeMemory = true;
                    pg.start();
                    while (parser.hasNextDataLine()) {
                        int i = 0;
                        for (final int globalFieldIndex : selectedInfoFieldGlobalIndices) {

                            final String stringFieldValue = parser.getStringFieldValue(globalFieldIndex);
                            stream.write(stringFieldValue);
                            if (i++ < selectedInfoFieldGlobalIndices.size()) {
                                stream.write("\t");
                            }
                        }

                        stream.write("\n");
                        parser.next();
                    }
                } catch (VCFParser.SyntaxException e) {
                    System.err.println("A syntax error was encountered when parsing input file. Details may be provided below. " + filename);
                    e.printStackTrace();
                }
            }
        } finally {
            if (outputFilename != null) {
                IOUtils.closeQuietly(stream);
            }
        }
    }

    private ColumnField findColumnField(String selectedFieldName, ColumnInfo... columns) {
        for (ColumnInfo col : columns) {
            ColumnField v = col.fields.find(selectedFieldName);
            if (v != null) {
                return v;
            }
        }
        return null;
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
        new VcfToTabMode().configure(args).execute();
    }
}
