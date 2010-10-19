/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.reads.ReadSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

/**
 * Concatenate compact reads files, count the number of reads, and
 * track the min and max sequence length of all of the reads.
 *
 * @author Kevin Dorff
 */
public class ConcatenateCompactReadsMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(ConcatenateCompactReadsMode.class);

    /**
     * The input files.
     */
    private List<File> inputFiles;

    /**
     * The output filename.
     */
    private String outputFilename;

    /**
     * sequences per chunck in the written file.
     */
    private int sequencePerChunk = 10000;

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "concatenate-compact-reads";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Concatenate compact reads files, count the "
            + "number of reads, and track the min and max sequence length of all of the reads.";

    private long numberOfReads;
    private int minReadLength = Integer.MAX_VALUE;
    private int maxReadLength = Integer.MIN_VALUE;
    private String optionalFilterExtension;

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
     * @throws IOException   error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        setInputFilenames(jsapResult.getStringArray("input"));
        outputFilename = jsapResult.getString("output");
        optionalFilterExtension = jsapResult.getString("optional-filter-extension");
        sequencePerChunk = jsapResult.getInt("sequence-per-chunk");
        return this;
    }

    /**
     * Actually perform the split of the compact reads file between
     * start and end position to the new compact reads file.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {
        if (inputFiles == null || inputFiles.size() == 0) {
            throw new IOException("--input not specified");
        }
        if (StringUtils.isBlank(outputFilename)) {
            throw new IOException("--output not specified");
        }

        final ReadsWriter writer = new ReadsWriter(new FileOutputStream(outputFilename));
        writer.setNumEntriesPerChunk(sequencePerChunk);
        final MutableString sequence = new MutableString();

        ReadsReader readsReader = null;
        numberOfReads = 0;
        minReadLength = Integer.MAX_VALUE;
        maxReadLength = Integer.MIN_VALUE;
        int removedByFilterCount = 0;
        try {
            for (final File inputFile : inputFiles) {
                readsReader = new ReadsReader(inputFile);
                String basename = FilenameUtils.removeExtension(inputFile.getPath());
                String filterFilename = FilenameUtils.concat(basename, optionalFilterExtension);
                File filterFile = new File(filterFilename);
                ReadSet readIndexFilter = null;
                if (filterFile.exists() && filterFile.canRead()) {
                    readIndexFilter = new ReadSet();
                    readIndexFilter.load(filterFile);
                    LOG.info(String.format("Loaded optional filter %s with %d elements. ",
                            filterFile, readIndexFilter.size()));
                }
                for (final Reads.ReadEntry readEntry : readsReader) {
                    // only concatenate if (1) there is no filter or (2) the read index is in the filter.
                    if (readIndexFilter == null || readIndexFilter.contains(readEntry.getReadIndex())) {

                        numberOfReads++;
                        if (readEntry.hasDescription()) {
                            writer.setDescription(readEntry.getDescription());
                        }
                        if (readEntry.hasReadIdentifier()) {
                            writer.setIdentifier(readEntry.getReadIdentifier());
                        }
                        if (readEntry.hasQualityScores()) {
                            final byte[] scores = ReadsReader.decodeQualityScores(readEntry);
                            if (scores != null) {
                                writer.setQualityScores(scores);
                            }
                        }
                        minReadLength = Math.min(minReadLength, readEntry.getReadLength());
                        maxReadLength = Math.max(maxReadLength, readEntry.getReadLength());
                        if (readEntry.hasSequence()) {
                            ReadsReader.decodeSequence(readEntry, sequence);
                            writer.setSequence(sequence);
                        } else {
                            writer.setSequence("");
                        }
                        writer.appendEntry();
                    } else {
                        removedByFilterCount++;
                    }
                }
                readsReader.close();
                readsReader = null;

            }
        } finally {
            writer.printStats(System.out);
            System.out.println("Number of reads=" + numberOfReads);
            System.out.println("Minimum Read Length=" + minReadLength);
            System.out.println("Maximum Read Length=" + maxReadLength);
            System.out.println("Reads removed by filter=" + removedByFilterCount);
            writer.close();
            if (readsReader != null) {
                readsReader.close();
            }
        }
    }

    /**
     * Add an input file.
     *
     * @param inputFile the input file to add.
     */
    public synchronized void addInputFile(final File inputFile) {
        if (inputFiles == null) {
            inputFiles = new LinkedList<File>();
        }
        this.inputFiles.add(inputFile);
    }

    /**
     * Clear the input files list.
     */
    public synchronized void clearInputFiles() {
        if (inputFiles != null) {
            inputFiles.clear();
        }
    }

    /**
     * Set the input filenames.
     *
     * @param inputFilenames the input filename
     */
    public synchronized void setInputFilenames(final String[] inputFilenames) {
        clearInputFiles();
        for (final String inputFilname : inputFilenames) {
            addInputFile(new File(inputFilname));
        }
    }

    /**
     * Get the input filenames.
     *
     * @return the input filenames
     */
    public synchronized String[] getInputFilenames() {
        if (inputFiles == null) {
            return new String[0];
        }
        final String[] array = new String[inputFiles.size()];
        int i = 0;
        for (final File inputFile : inputFiles) {
            array[i++] = inputFile.toString();
        }
        return array;
    }

    /**
     * Set the output filename.
     *
     * @param outputFilename the output filename
     */
    public void setOutputFilename(final String outputFilename) {
        this.outputFilename = outputFilename;
    }

    /**
     * Get the output filename.
     *
     * @return the output filename
     */
    public String getOutputFilename() {
        return this.outputFilename;
    }

    /**
     * The number of reads after the concatenate.
     *
     * @return number of reads
     */
    public long getNumberOfReads() {
        return numberOfReads;
    }

    /**
     * The minimum sequence length of all of the reads.
     *
     * @return minimum sequence length
     */
    public int getMinReadLength() {
        return minReadLength;
    }

    /**
     * The maximum sequence length of all of the reads.
     *
     * @return maximum sequence length
     */
    public int getMaxReadLength() {
        return maxReadLength;
    }

    /**
     * Main mode for splitting compact reads files from a start position
     * to and end position.
     *
     * @param args command line arguments
     * @throws java.io.IOException IO error
     * @throws com.martiansoftware.jsap.JSAPException
     *                             command line parsing error.
     */
    public static void main(final String[] args) throws IOException, JSAPException {
        new ConcatenateCompactReadsMode().configure(args).execute();
    }
}
