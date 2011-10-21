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
import edu.cornell.med.icb.io.TsvToFromMap;
import edu.cornell.med.icb.iterators.TsvLineIterator;
import edu.cornell.med.icb.maps.LinkedHashToMultiTypeMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectLinkedOpenHashMap;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Read the data from TSV files to determine the the column types (double/integer/string).
 * Write a .colinfo file detailing the column names and types.
 */
public class TabToColumnInfoMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TabToColumnInfoMode.class);

    /** The possible column types. */
    public enum ColumnTypes {
        STRING,
        DOUBLE,
        INTEGER,
        UNKNOWN
    }

    private static final Pattern DOUBLE_PATTERN =
            Pattern.compile("[\\-\\+]?\\d+(\\.\\d*)?(e-?\\d+)?|nan|-infinity|\\+infinity|infinity|inf|-inf|\\+inf");
    private static final Pattern INTEGER_PATTERN =
            Pattern.compile("-?\\d+?");
    private static final Set<String> SPECIAL_DOUBLES = new HashSet<String>();
    static {
        // Only use lowercase for these.
        SPECIAL_DOUBLES.add("nan");
        SPECIAL_DOUBLES.add("-infinity");
        SPECIAL_DOUBLES.add("+infinity");
        SPECIAL_DOUBLES.add("infinity");
        SPECIAL_DOUBLES.add("-inf");
        SPECIAL_DOUBLES.add("+inf");
        SPECIAL_DOUBLES.add("inf");
    }

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "tab-to-column-info";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Read the data from TSV files to determine the the column types " +
            "(double/integer/string). Write a .colinfo file detailing the column names and types.";

    /**
     * The input files.
     */
    private List<File> inputFiles;

    /**
     * The number of input lines to read or set to a value less than or equal to 0 to read the entire file.
     */
    private int numberOfLinesToProcess = 10000;

    /**
     * Set to true if you don't want to output any data, just process the file and leave
     * the results in data structures in the object, good for API calls to this class
     * or testing.
     */
    private boolean noOutput;

    /**
     * If we are running in API mode. In API mode exceptions will be thrown instead of System.exit() when there
     * are problems. Running configure() to parse command line options automatically turns apiMode off.
     */
    private boolean apiMode = true;

    /**
     * The filename to a (map of the columnName to the column type).
     * The order of the filename keys is preserved.
     */
    private final Map<String, Map<String, ColumnTypes>> filenameToDetailsMap =
            new Object2ObjectLinkedOpenHashMap<String, Map<String, ColumnTypes>>();

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
        apiMode = false;
        return this;
    }

    /**
     * Add an input file.
     *
     * @param inputFile the input file to add.
     */
    public void addInputFile(final File inputFile) {
        if (inputFiles == null) {
            inputFiles = new LinkedList<File>();
        }
        inputFiles.add(inputFile);
    }

    /**
     * Clear the input files list.
     */
    public void clearInputFiles() {
        if (inputFiles != null) {
            inputFiles.clear();
        }
    }

    /**
     * Set the input filenames.
     *
     * @param inputFilenames the input filename
     */
    public void setInputFilenames(final String[] inputFilenames) {
        clearInputFiles();
        for (final String inputFilename : inputFilenames) {
            addInputFile(new File(inputFilename));
        }
    }

    /**
     * Get if noOutput mode is enabled. If true, no output files are written. Default is false.
     * @return the value of noOutput
     */
    public boolean isNoOutput() {
        return noOutput;
    }

    /**
     * Set if noOutput mode is enabled. If true, no output files are written. Default is false.
     * @param noOutput the new value of noOutput
     */
    public void setNoOutput(final boolean noOutput) {
        this.noOutput = noOutput;
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

    public Map<String, Map<String, ColumnTypes>> getFilenameToDetailsMap() {
        return filenameToDetailsMap;
    }

    /**
     * Gather column info for provided input files.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {
        if (inputFiles == null || inputFiles.isEmpty()) {
            if (apiMode) {
                throw new IOException("No input files provided");
            } else {
                System.err.println("No input files provided");
                System.exit(1);
            }
        }
        try {
            for (final File inputFile : inputFiles) {
                processOneFile(inputFile);
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
     * @param inputFile the file to process
     * @throws IOException an error processing the file.
     */
    public void processOneFile(final File inputFile) throws IOException {
        final String filename = inputFile.toString();
        if (!inputFile.exists()) {
            throw new IOException("File not found " + filename);
        }

        TsvLineIterator in = null;
        int numProcessed = 0;
        PrintWriter out = null;
        final boolean quitEarly = numberOfLinesToProcess > 0;
        try {

            if (filenameToDetailsMap.get(filename) != null) {
                // Don't process the same file twice.
                return;
            }

            final Map<String, ColumnTypes> columnToType = new Object2ObjectLinkedOpenHashMap<String, ColumnTypes>();
            filenameToDetailsMap.put(filename, columnToType);

            // Get the column headers, initialize columnToType
            final TsvToFromMap tsvDetails = TsvToFromMap.createFromTsvFile(inputFile);
            for (final String columnName : tsvDetails.getColumnHeaders()) {
                columnToType.put(columnName, ColumnTypes.UNKNOWN);
            }

            // Iterate over the file.
            in = new TsvLineIterator(inputFile, tsvDetails);
            for (final LinkedHashToMultiTypeMap<String> lineMap : in) {
                for (final Map.Entry<String, String> entry : lineMap.entrySet()) {
                    final String columnName = entry.getKey();
                    final String columnValue = entry.getValue();
                    final ColumnTypes prevColumnType = columnToType.get(columnName);
                    if (prevColumnType != ColumnTypes.STRING) {
                        // If this column was previously known to be STRING, nothing to do
                        // But since it's not, we need to observe the value.
                        final ColumnTypes columnType = typeFromValue(columnValue);
                        if (columnType == ColumnTypes.UNKNOWN) {
                            // Unknown means the value was empty. This shouldn't have ANY
                            // effect on the column type.
                        } else if (columnType == ColumnTypes.STRING) {
                            // If we see any strings in a column, the column is forced to be a String
                            // without respect to any other values
                            columnToType.put(columnName, ColumnTypes.STRING);
                        } else if (prevColumnType == ColumnTypes.UNKNOWN) {
                            // Uninitialized, set to the new type
                            columnToType.put(columnName, columnType);
                        } else {
                            if (columnType == ColumnTypes.DOUBLE && prevColumnType == ColumnTypes.INTEGER) {
                                // We need to promote this integer to a double
                                columnToType.put(columnName, ColumnTypes.DOUBLE);
                            }
                        }
                    }
                }

                if (quitEarly && ++numProcessed == numberOfLinesToProcess) {
                    // Quit early if requested
                    break;
                }
            }

            if (!noOutput) {
                out = new PrintWriter(filename + ".colinfo");
                int columnNumber = 0;
                for (final Map.Entry<String, ColumnTypes> entry : columnToType.entrySet()) {
                    if (columnNumber == 0) {
                        out.println("column-name\tcolumn-type");
                    }
                    out.print(entry.getKey());
                    out.print('\t');
                    if (entry.getValue() == ColumnTypes.UNKNOWN) {
                        // If there are no values in the column, we'll just assume the column to be a string
                        out.print(ColumnTypes.STRING.toString().toLowerCase());
                    } else {
                        out.print(entry.getValue().toString().toLowerCase());
                    }

                    out.print(IOUtils.LINE_SEPARATOR);
                    columnNumber++;
                }
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (out != null) {
                IOUtils.closeQuietly(out);
            }
        }
    }


    /**
     * Check the type of value, is it an integer, double, or a string.
     * @param value the value to check
     * @return the column type that is suitable for value
     */
    public static ColumnTypes typeFromValue(final String value) {
        if (value == null) {
            return ColumnTypes.UNKNOWN;
        }
        final String v = value.toLowerCase().trim();
        if (v.isEmpty()) {
            return ColumnTypes.UNKNOWN;
        }
        if (INTEGER_PATTERN.matcher(v).matches()) {
            // We first test for Integer because double would match it but not vice versa.
            try {
                // We SHOULD have a valid Integer value, but let's make sure.
                final Integer i = Integer.valueOf(v);
                return ColumnTypes.INTEGER;
            } catch (java.lang.NumberFormatException e) {
                // It isn't a double. Move on to the next case.
                System.err.println("Converting " + v + " to integer caused an exception");
            }
        }
        if (DOUBLE_PATTERN.matcher(v).matches()) {
            try {
                // Check for "special" values.
                if (SPECIAL_DOUBLES.contains(v)) {
                    return ColumnTypes.DOUBLE;
                }
                // We SHOULD have a valid Double, but let's make sure.
                final Double d = Double.valueOf(v);
                return ColumnTypes.DOUBLE;
            } catch (java.lang.NumberFormatException e) {
                // It isn't a double. Move on to the next case.
                System.err.println("Converting " + v + " to double caused an exception");
            }
        }
        // Wasn't an integer or double, it must be a string.
        return ColumnTypes.STRING;
    }
}
