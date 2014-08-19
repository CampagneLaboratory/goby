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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.io.TsvToFromMap;
import edu.cornell.med.icb.iterators.TsvLineIterator;
import edu.cornell.med.icb.maps.LinkedHashToMultiTypeMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectLinkedOpenHashMap;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Read the data from TSV files to determine the the column types (Float/Integer/String).
 * Write a .colinfo file detailing the column names and types.
 */
public class TabToColumnInfoMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TabToColumnInfoMode.class);

    private static final Pattern FLOAT_PATTERN =
            Pattern.compile("[\\-\\+]?\\d+(\\.\\d*)?(e-?\\d+)?|nan|-infinity|\\+infinity|infinity|inf|-inf|\\+inf");
    private static final Pattern INTEGER_PATTERN =
            Pattern.compile("-?\\d+?");
    private static final Set<String> SPECIAL_FLOATS = new HashSet<String>();
    static {
        // Only use lowercase for these.
        SPECIAL_FLOATS.add("nan");
        SPECIAL_FLOATS.add("-infinity");
        SPECIAL_FLOATS.add("+infinity");
        SPECIAL_FLOATS.add("infinity");
        SPECIAL_FLOATS.add("-inf");
        SPECIAL_FLOATS.add("+inf");
        SPECIAL_FLOATS.add("inf");
    }

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "tab-to-column-info";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Read the data from TSV files to determine the the column types " +
            "(Float/Integer/String). Write a .colinfo file detailing the column names and types.";

    /**
     * The input filenames.
     */
    private Set<String> inputFilenames;

    /**
     * The number of input lines to read or set to a value less than or equal to 0 to read the entire file.
     */
    private int numberOfLinesToProcess = 10000;

    /**
     * Set to true if you don't want to output any data, just process the file and leave
     * the results in data structures in the object, good for API calls to this class
     * or testing.
     */
    private boolean createCache;

    /**
     * If we are running in API mode. In API mode exceptions will be thrown instead of System.exit() when there
     * are problems. Running configure() to parse command line options automatically turns apiMode off.
     */
    private boolean apiMode = true;

    /**
     * If the .colinfo file already exists, if this is true it will be used. By default, from the command line
     * it will ALWAYS re-create the file. This is for API use.
     */
    private boolean readFromCache;

    /**
     * Output the results to stdout.
     */
    private boolean display;

    /**
     * Verbose.
     */
    private boolean verbose;

    /**
     * The filename to a (map of the columnName to the column type).
     * The order of the filename keys is preserved.
     */
    private final Map<String, Map<String, ColumnType>> filenameToDetailsMap =
            new Object2ObjectLinkedOpenHashMap<String, Map<String, ColumnType>>();

    /**
     * Mode name.
     * @return Mode name.
     */
    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    /**
     * Mode description.
     * @return Mode description.
     */
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
        setInputFilenames(jsapResult.getStringArray("input"));
        numberOfLinesToProcess = jsapResult.getInt("number-of-lines");
        createCache = !jsapResult.getBoolean("do-not-create-cache");
        display = jsapResult.getBoolean("display");
        readFromCache = jsapResult.getBoolean("read-from-cache");
        verbose = jsapResult.getBoolean("verbose");
        apiMode = false;
        return this;
    }

    /**
     * Add an input filename.
     *
     * @param inputFilename the input filename to add.
     */
    public void addInputFilename(final String inputFilename) {
        if (inputFilenames == null) {
            inputFilenames = new LinkedHashSet<String>();
        }
        inputFilenames.add(inputFilename);
    }

    /**
     * Add an input file.
     *
     * @param inputFile the input file to add.
     */
    public void addInputFile(final File inputFile) {
        addInputFilename(inputFile.toString());
    }

    /**
     * Clear the input files list.
     */
    public void clearInputFilenames() {
        if (inputFilenames != null) {
            inputFilenames.clear();
        }
    }

    /**
     * Set the input filenames.
     *
     * @param inputFilenames the input filename
     */
    public void setInputFilenames(final String[] inputFilenames) {
        clearInputFilenames();
        for (final String inputFilename : inputFilenames) {
            addInputFilename(inputFilename);
        }
    }

    /**
     * Set the input filenames.
     *
     * @param inputFiles the input filename
     */
    public void setInputFiles(final File[] inputFiles) {
        clearInputFilenames();
        for (final File inputFile : inputFiles) {
            addInputFile(inputFile);
        }
    }

    /**
     * Get if createCache mode is enabled. If true, no output files are written. Default is false.
     * @return the value of createCache
     */
    public boolean isCreateCache() {
        return createCache;
    }

    /**
     * Set if createCache mode is enabled. If true, no output files are written. Default is false.
     * @param createCache the new value of createCache
     */
    public void setCreateCache(final boolean createCache) {
        this.createCache = createCache;
    }

    /**
     * Output the results to stdout.
     * @return if display enabled
     */
    public boolean isDisplay() {
        return display;
    }

    /**
     * Output the results to stdout.
     * @param display new value for display
     */
    public void setDisplay(final boolean display) {
        this.display = display;
    }

    /**
     * Get the number of input lines to read or set to <= 0 to read the entire file.
     * @return the value of numberOfLinesToProcess
     */
    public int getNumberOfLinesToProcess() {
        return numberOfLinesToProcess;
    }

    /**
     * Set the number of input lines to read or set to <= 0 to read the entire file.
     * @param numberOfLinesToProcess the new value of numberOfLinesToProcess
     */
    public void setNumberOfLinesToProcess(final int numberOfLinesToProcess) {
        this.numberOfLinesToProcess = numberOfLinesToProcess;
    }

    /**
     * If true, the [filename].colinfo file already exists, if this is true it will be used.
     * @return if readFromCache is enabled.
     */
    public boolean isReadFromCache() {
        return readFromCache;
    }

    /**
     * If true, the [filename].colinfo file already exists, if this is true it will be used.
     * @param readFromCache if readFromCache is enabled.
     */
    public void setReadFromCache(final boolean readFromCache) {
        this.readFromCache = readFromCache;
    }

    /**
     * Verbose.
     * @return verbose value
     */
    public boolean isVerbose() {
        return verbose;
    }

    /**
     * Verbose.
     * @param verbose new verbose value
     */
    public void setVerbose(final boolean verbose) {
        this.verbose = verbose;
    }

    /**
     * Get the map of filenames -> (columnName -> columnType map).
     * @return the map of files to column details.
     */
    public Map<String, Map<String, ColumnType>> getFilenameToDetailsMap() {
        return filenameToDetailsMap;
    }

    /**
     * Get a specific column details based on the order of the input filenames given.
     * @param index the index of results to retrieve
     * @return the details for the specified index
     * @throws IndexOutOfBoundsException if index provided is too large
     */
    public Map<String, ColumnType> getDetailsAtIndex(final int index) throws IndexOutOfBoundsException {
        if (index < 0) {
            throw new IndexOutOfBoundsException("Index must be >= 0");
        }
        final int maxIndex = filenameToDetailsMap.size() - 1;
        if (index > maxIndex) {
            throw new IndexOutOfBoundsException("Index" + index + " should have been <= " + maxIndex);
        }

        final String[] inputFilenamesAsArray = inputFilenames.toArray(new String[inputFilenames.size()]);
        return filenameToDetailsMap.get(inputFilenamesAsArray[index]);
    }

    /**
     * Gather column info for provided input files.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {
        if (inputFilenames == null || inputFilenames.isEmpty()) {
            if (apiMode) {
                throw new IOException("No input files provided");
            } else {
                System.err.println("No input files provided");
                System.exit(1);
            }
        }
        try {
            for (final String inputFilename : inputFilenames) {
                processOneFile(inputFilename);
                if (display) {
                    displayToStdout(inputFilename, filenameToDetailsMap.get(inputFilename));
                }
            }
        } catch (IOException e) {
            LOG.error("Error processing", e);
            if (apiMode) {
                throw e;
            } else {
                System.exit(2);
            }
        }
    }

    /**
     * Process one file
     * @param filename the filename to process
     * @throws IOException an error processing the file.
     */
    public void processOneFile(final String filename) throws IOException {
        final File inputFile = new File(filename);
        if (!inputFile.exists()) {
            throw new IOException("File not found " + filename);
        }

        if (readFromCache) {
            if (readFromCache(filename)) {
                return;
            }
        }

        TsvLineIterator in = null;
        int numProcessed = 0;
        final boolean quitEarly = numberOfLinesToProcess > 0;
        try {

            if (filenameToDetailsMap.get(filename) != null) {
                // Don't process the same file twice.
                return;
            }

            final Map<String, ColumnType> columnToType = new Object2ObjectLinkedOpenHashMap<String, ColumnType>();
            filenameToDetailsMap.put(filename, columnToType);

            // Get the column headers, initialize columnToType
            final TsvToFromMap tsvDetails = TsvToFromMap.createFromTsvFile(inputFile);
            tsvDetails.setLenientColumnCount(true);
            for (final String columnName : tsvDetails.getColumnHeaders()) {
                columnToType.put(columnName, ColumnType.Unknown);
            }

            // Iterate over the file.
            in = new TsvLineIterator(inputFile, tsvDetails);

            for (final LinkedHashToMultiTypeMap<String> lineMap : in) {
                for (final Map.Entry<String, String> entry : lineMap.entrySet()) {
                    final String columnName = entry.getKey();
                    final String columnValue = entry.getValue();
                    final ColumnType prevColumnType = columnToType.get(columnName);
                    if (prevColumnType != ColumnType.String) {
                        // If this column was previously known to be STRING, nothing to do
                        // But since it's not, we need to observe the value.
                        final ColumnType columnType = typeFromValue(columnValue);
                        if (columnType == ColumnType.Unknown) {
                            // Unknown means the value was empty. This shouldn't have ANY
                            // effect on the column type.
                        } else if (columnType == ColumnType.String) {
                            // If we see any strings in a column, the column is forced to be a String
                            // without respect to any other values
                            columnToType.put(columnName, ColumnType.String);
                        } else if (prevColumnType == ColumnType.Unknown) {
                            // Uninitialized, set to the new type
                            columnToType.put(columnName, columnType);
                        } else {
                            if (columnType == ColumnType.Float && prevColumnType == ColumnType.Integer) {
                                // We need to promote this integer to a float
                                columnToType.put(columnName, ColumnType.Float);
                            }
                        }
                    }
                }

                if (quitEarly && ++numProcessed == numberOfLinesToProcess) {
                    // Quit early if requested
                    break;
                }
            }
            if (createCache) {
                output(filename + ".colinfo", columnToType);
            }
        } finally {
            if (in != null) {
                in.close();
            }
        }
    }

    private void output(final String cacheFilename, final Map<String, ColumnType> columnToType)
            throws FileNotFoundException {
        PrintWriter out = null;
        try {
            out = new PrintWriter(cacheFilename);
            int columnNumber = 0;
            for (final Map.Entry<String, ColumnType> entry : columnToType.entrySet()) {
                if (columnNumber == 0) {
                    out.println("id\tcolumn-name\tcolumn-type");
                }
                out.print("col_" + columnNumber);
                out.print('\t');
                out.print(entry.getKey());
                out.print('\t');
                if (entry.getValue() == ColumnType.Unknown) {
                    // If there are no values in the column, we'll just assume the column to be a string
                    out.print(ColumnType.String.toString());
                } else {
                    out.print(entry.getValue().toString());
                }

                out.print(IOUtils.LINE_SEPARATOR);
                columnNumber++;
            }
        } finally {
            if (out != null) {
                IOUtils.closeQuietly(out);
            }
        }
    }

    private void displayToStdout(final String filename, final Map<String, ColumnType> columnToType) {
        if (display) {
            System.out.println("filename=" + filename);
            for (final Map.Entry<String, ColumnType> entry : columnToType.entrySet()) {
                System.out.println("  " + entry.getKey() + "=" + entry.getValue());
            }
        }
    }

    private boolean readFromCache(final String inputTsvFilename) throws IOException {
        final File inputTsvFile = new File(inputTsvFilename);
        final File inputColInfoFile = new File(inputTsvFile + ".colinfo");
        if (!inputColInfoFile.exists()) {
            if (verbose) {
                System.err.println("#Cached .colinfo file [" + inputColInfoFile.toString() + "] not found");
            }
            return false;
        }
        final TsvToFromMap columns = TsvToFromMap.createFromTsvFile(inputColInfoFile);
        final List<String> columnNames = columns.getColumnHeaders();
        if (columnNames.size() != 3) {
            if (verbose) {
                System.err.println("#Cached .colinfo file contains the wrong number of columns");
            }
            return false;
        }
        // Make sure .colinfo file contains the right columns
        if (!(columnNames.contains("id") &&
                columnNames.contains("column-name") &&
                columnNames.contains("column-type"))) {
            if (verbose) {
                System.err.println("#Cached .colinfo contains the wrong column names");
            }
            return false;
        }

        final Map<String, ColumnType> columnToType = new Object2ObjectLinkedOpenHashMap<String, ColumnType>();
        filenameToDetailsMap.put(inputTsvFilename, columnToType);

        // Read the column types into the map
        for (final LinkedHashToMultiTypeMap<String> lineMap : new TsvLineIterator(inputColInfoFile, columns)) {
            columnToType.put(lineMap.getString("column-name"),
                    ColumnType.valueOf(lineMap.getString("column-type")));
        }
        if (columnToType.isEmpty()) {
            // No columns found in .colinfo file.
            if (verbose) {
                System.err.println("#Cached .colinfo contains no data");
            }
            return false;
        }

        // Read the columns from the original TSV file
        final List<String> origTsvColumnsNames = TsvToFromMap.createFromTsvFile(inputTsvFile).getColumnHeaders();
        if (origTsvColumnsNames.size() != columnToType.size()) {
            if (verbose) {
                System.err.println("#Cached .colinfo contains the wrong number of data rows");
            }
            return false;
        }
        for (final String columnName : origTsvColumnsNames) {
            if (!columnToType.containsKey(columnName)) {
                // Column missing
                filenameToDetailsMap.remove(inputTsvFilename);
                if (verbose) {
                    System.err.println("#Cached .colinfo data rows contain the wrong column names");
                }
                return false;
            }
        }
        if (verbose) {
            System.out.println("#Column info read from cache.");
        }
        return true;
    }

    /**
     * Check the type of value, is it an Integer, Float, or a String.
     * @param value the value to check
     * @return the column type that is suitable for value
     */
    public ColumnType typeFromValue(final String value) {
        if (value == null) {
            return ColumnType.Unknown;
        }
        final String v = value.toLowerCase().trim();
        if (v.isEmpty()) {
            return ColumnType.Unknown;
        }
        if (INTEGER_PATTERN.matcher(v).matches()) {
            // We first test for Integer because Float would match it but not vice versa.
            try {
                // We SHOULD have a valid Integer value, but let's make sure.
                final Integer i = Integer.valueOf(v);
                return ColumnType.Integer;
            } catch (java.lang.NumberFormatException e) {
                // It isn't a Float. Move on to the next case.
                if (verbose) {
                    System.err.println("#Converting " + v + " to integer caused an exception");
                }
            }
        }
        if (FLOAT_PATTERN.matcher(v).matches()) {
            try {
                // Check for "special" values.
                if (SPECIAL_FLOATS.contains(v)) {
                    return ColumnType.Float;
                }
                // We SHOULD have a valid Double, but let's make sure.
                final Double d = Double.valueOf(v);
                return ColumnType.Float;
            } catch (java.lang.NumberFormatException e) {
                // It isn't a double. Move on to the next case.
                if (verbose) {
                    System.err.println("#Converting " + v + " to double caused an exception");
                }
            }
        }
        // Wasn't an Integer or Float, it must be a String.
        return ColumnType.String;
    }
}
