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
import edu.cornell.med.icb.goby.readers.vcf.ColumnField;
import edu.cornell.med.icb.goby.readers.vcf.ColumnInfo;
import edu.cornell.med.icb.goby.readers.vcf.Columns;
import edu.cornell.med.icb.goby.readers.vcf.VCFParser;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import edu.cornell.med.icb.goby.util.DoInParallel;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Replacement for vcf-subset tool. This mode provides a replacement for the vcf-tools vcf-subset utility.
 * In contrast to vcf-subset (vcf-tools), this implementation can extract specific samples from GB of data
 * from the 1000g project in a few minutes (where vcf-subset would require hours of compute time for each
 * chromosome file).
 *
 * @author Fabien Campagne
 */
public class VCFSubsetMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(VCFSubsetMode.class);

    /**
     * The output filename.
     */
    private String outputFilename;


    /**
     * The mode name.
     */
    private static final String MODE_NAME = "vcf-subset";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Extract specific samples from large VCF files.";

    private int[] chromosomeFieldIndex;
    private int[] positionFieldIndex;

    private ObjectArraySet<String> sampleIdsSelected;
    private int[] sampleIndexToDestinationIndex;
    private String[] inputFilenames;
    private boolean doInParallel;

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
        inputFilenames = jsapResult.getStringArray("input");
        outputFilename = jsapResult.getString("output");
        String[] columns = jsapResult.getStringArray("column");
        this.sampleIdsSelected = new ObjectArraySet<String>(columns);
        doInParallel=jsapResult.getBoolean("parallel");
        if (sampleIdsSelected.size() == 0) {
            System.err.println("You must select at least one column.");
            System.exit(1);
        }
        return this;
    }


    /**
     * Compare VCF files.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {
        if (inputFilenames == null || inputFilenames.length == 0) {
            throw new IOException("--input not specified");
        }
        if (StringUtils.isBlank(outputFilename)) {
            throw new IOException("--output not specified");
        }
        final int numInputFiles = inputFilenames.length;

        int parserIndex = 0;
        chromosomeFieldIndex = new int[numInputFiles];
        positionFieldIndex = new int[numInputFiles];
        DoInParallel loop = new DoInParallel() {
            @Override
            public void action(DoInParallel forDataAccess, String inputBasename, int loopIndex) {
                try {
                    processOneFile(new File(inputBasename));
                } catch (IOException e) {
                    e.printStackTrace(System.err);
                    System.exit(1);
                }
            }
        };
        try {
            loop.execute(doInParallel, inputFilenames);
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(1);
        }
    }

    private void processOneFile(File inputFile) throws IOException {
        System.out.printf("Preparing to process %s..%n", inputFile);
        final VCFParser parser = new VCFParser(inputFile.getPath());
        final Columns columns = new Columns();
        final ObjectArrayList<String> sampleIdList = new ObjectArrayList<String>();
        boolean[] includeField = null;
        try {
            parser.readHeader();
            final Columns fileColumns = parser.getColumns();

            includeField = new boolean[parser.countAllFields()];

            for (final ColumnInfo col : fileColumns) {
                if (columns.find(col.getColumnName()) == null) {
                    columns.add(col);
                    if (col.useFormat) {
                        final String sampleId = col.getColumnName();
                        if (!sampleIdList.contains(sampleId) && includeSampleId(sampleId)) {
                            for (final ColumnField field : col.fields) {
                                includeField[field.globalFieldIndex] = true;
                            }
                            sampleIdList.add(sampleId);

                        }
                    } else {
                        for (final ColumnField field : col.fields) {
                            includeField[field.globalFieldIndex] = true;
                        }
                    }
                }
            }
        } catch (VCFParser.SyntaxException e) {
            e.printStackTrace();
            System.exit(1);
        }
        final String inputFilename = removeExtensions(inputFile);

        final int chromosomeFieldIndex = getGlobalFieldIndex(columns, "CHROM");
        final int positionFieldIndex = getGlobalFieldIndex(columns, "POS");
        final int idFieldIndex = getGlobalFieldIndex(columns, "ID");
        final int refFieldIndex = getGlobalFieldIndex(columns, "REF");
        final int altFieldIndex = getGlobalFieldIndex(columns, "ALT");
        final int qualFieldIndex = getGlobalFieldIndex(columns, "QUAL");
        final int filterFieldIndex = getGlobalFieldIndex(columns, "FILTER");

        final IntSet infoFieldGlobalIndices = new IntArraySet();
        sampleIndexToDestinationIndex = new int[parser.countAllFields()];
        for (final ColumnField infoField : parser.getColumns().find("INFO").fields) {
            infoFieldGlobalIndices.add(infoField.globalFieldIndex);
        }

        final IntSet formatFieldGlobalIndices = new IntArraySet();
        final Int2IntMap globalIndexToSampleIndex = new Int2IntOpenHashMap();
        int sampleIndex = 0;

        for (ColumnInfo col : columns) {
            if (col.useFormat) {

                for (ColumnField field : col.fields) {
                    String sampleId = col.getColumnName();
                    if (sampleIdList.contains(sampleId)) {
                        sampleIndexToDestinationIndex[sampleIndex] = sampleIdList.indexOf(sampleId);
                    }
                    globalIndexToSampleIndex.put(field.globalFieldIndex, sampleIndex++);
                    formatFieldGlobalIndices.add(field.globalFieldIndex);
                }
            }
        }

        int infoFieldIndex = 0;
        int formatFieldCount = 0;
        int previousSampleIndex = -1;

        // transfer the reduced schema to the output writer:
        VCFWriter writer = new VCFWriter(new FileWriter(inputFilename + outputFilename + ".vcf"));
        writer.defineSchema(columns);
        writer.defineSamples(sampleIdList.toArray(new String[sampleIdList.size()]));
        writer.writeHeader();
        System.out.printf("Loading %s..%n", inputFilename);
        int index = 0;
        final ProgressLogger pg = new ProgressLogger(LOG);
        pg.displayFreeMemory = true;
        pg.start();
        final int fieldCount = parser.countAllFields();
        final IntSet fieldsToTraverse = new IntArraySet();
        fieldsToTraverse.addAll(infoFieldGlobalIndices);
        fieldsToTraverse.addAll(formatFieldGlobalIndices);
        for (int i = 0; i < filterFieldIndex; i++) {
            fieldsToTraverse.add(i);
        }
        while (parser.hasNextDataLine()) {
            final String format = parser.getStringColumnValue(columns.find("FORMAT").columnIndex);
            final String[] formatTokens = format.split(":");
            final int numFormatFields = columns.find("FORMAT").fields.size();

            int formatFieldIndex = 0;
            infoFieldIndex = 0;
            for (final int globalFieldIndex : fieldsToTraverse) {
                final String value = parser.getStringFieldValue(globalFieldIndex);
                if (globalFieldIndex == chromosomeFieldIndex) {
                    writer.setChromosome(value);
                } else if (globalFieldIndex == positionFieldIndex) {
                    writer.setPosition(Integer.parseInt(value));
                } else if (globalFieldIndex == idFieldIndex) {
                    writer.setId(value);
                } else if (globalFieldIndex == refFieldIndex) {
                    writer.setReferenceAllele(value);
                } else if (globalFieldIndex == altFieldIndex) {
                    writer.setAlternateAllele(value);
                } else if (globalFieldIndex == qualFieldIndex) {
                    writer.setQual(value);
                } else if (globalFieldIndex == filterFieldIndex) {
                    writer.setFilter(value);
                }
                if (infoFieldGlobalIndices.contains(globalFieldIndex)) {
                    writer.setInfo(infoFieldIndex++, value);
                }

                if (formatFieldGlobalIndices.contains(globalFieldIndex)) {

                    if (formatFieldIndex < formatTokens.length) {
                        if (!"".equals(formatTokens[formatFieldIndex])) {
                            if (includeField[globalFieldIndex]) {
                                final int destinationSampleIndex = sampleIndexToDestinationIndex[globalIndexToSampleIndex.get(globalFieldIndex)];
                                writer.setSampleValue(formatTokens[formatFieldIndex], destinationSampleIndex, value);
                            }
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
            }
            parser.next();
            pg.lightUpdate();
            writer.writeRecord();
        }
        pg.stop("Done with file " + inputFilename);
        parser.close();

    }

    /**
     * Get the VALUE global field index for a column. If no VALUE field exist, but the column has
     * a single field, return that field. Fail (returns -1) if the column has more than one field.
     *
     * @param columns
     * @param columnName
     * @return
     */
    private int getGlobalFieldIndex(final Columns columns, final String columnName) {
        final ColumnInfo col = columns.find(columnName);
        final ColumnField value = col.getField("VALUE");
        if (value != null) {
            return value.globalFieldIndex;
        } else {
            if (col.fields.size() == 1) {
                return col.fields.get(0).globalFieldIndex;
            } else {
                System.err.println("Unable to obtain default field name for column " + columnName);
                return -1;
            }
        }
    }

    private boolean includeSampleId(String sampleId) {
        return sampleIdsSelected.contains(sampleId);
    }

    int fileIndex = 1;

    private String removeExtensions(File inputFile) {
        String filename = inputFile.getName();
        if (filename.endsWith(".vcf.gz")) {
            return FilenameUtils.removeExtension(FilenameUtils.removeExtension(inputFile.getName()));
        } else if (filename.endsWith(".vcf")) {
            return FilenameUtils.removeExtension(inputFile.getName());
        } else {
            return "output" + (fileIndex++);
        }
    }


    /**
     * @param args command line arguments
     * @throws java.io.IOException IO error
     * @throws com.martiansoftware.jsap.JSAPException
     *                             command line parsing error.
     */
    public static void main(final String[] args) throws IOException, JSAPException {
        new VCFSubsetMode().configure(args).execute();
    }

}