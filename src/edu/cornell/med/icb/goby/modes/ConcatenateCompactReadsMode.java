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
import edu.cornell.med.icb.goby.compression.MessageChunksWriter;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.reads.ReadSet;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.FileInputStream;
import java.util.LinkedList;
import java.util.List;
import java.nio.channels.FileChannel;

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

    /**
     * Display verbose output.
     */
    private boolean quickConcat;

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
        quickConcat = jsapResult.getBoolean("quick-concat", false);
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

        if (quickConcat) {
            performQuickConcat();
        } else {

            final ReadsWriter writer = new ReadsWriter(new FileOutputStream(outputFilename));
            writer.setNumEntriesPerChunk(sequencePerChunk);
            final MutableString sequence = new MutableString();

            ReadsReader readsReader = null;
            numberOfReads = 0;
            minReadLength = Integer.MAX_VALUE;
            maxReadLength = Integer.MIN_VALUE;
            int removedByFilterCount = 0;
            try {
                final ProgressLogger progress = new ProgressLogger();
                progress.start("concatenating files");
                progress.displayFreeMemory = true;
                progress.expectedUpdates = inputFiles.size();
                progress.start();
                for (final File inputFile : inputFiles) {

                    readsReader = new ReadsReader(inputFile);
                    String basename = FilenameUtils.removeExtension(inputFile.getPath());
                    String filterFilename = basename + optionalFilterExtension;
                    File filterFile = new File(filterFilename);
                    ReadSet readIndexFilter = null;
                    if (filterFile.exists() && filterFile.canRead()) {
                        readIndexFilter = new ReadSet();
                        readIndexFilter.load(filterFile);
                        LOG.info(String.format("Loaded optional filter %s with %d elements. ",
                                filterFile, readIndexFilter.size()));
                    } else{
                        if (optionalFilterExtension!=null) {
                            LOG.info("Could not locate filter for filename "+filterFilename);
                        }
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
                    progress.update();
                }
                progress.stop();
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
    }

    /**
     * This version does a quick concat. It does NO filtering. It gathers no stats,
     * but, will quickly concat multiple compact-reads files together using NIO.
     * It should be noted that this method is >MUCH< faster.
     * Copy all of the input files except the last MessageChunksWriter.DELIMITER_LENGTH
     * bytes of the first n-1 input files and the entire last input file
     * to the output file.
     * @throws IOException
     */
    private void performQuickConcat() throws IOException {
        System.out.println("quick concatenating files");
        File outputFile = new File(outputFilename);
        if (outputFile.exists()) {
            System.err.println("The output file already exists. Please delete it before running concat.");
            return;
        }
        outputFile.createNewFile();

        FileChannel input = null;
        FileChannel output = null;
        long maxChunkSize = 10*1024*1024; // 10 megabytes at a chunk
        try {
            output = new FileOutputStream(outputFile).getChannel();
            int lastFileNumToCopy = inputFiles.size() - 1;
            int curFileNum = 0;
            for (final File inputFile : inputFiles) {
                System.out.printf("Reading from %s%n", inputFile);
                input = new FileInputStream(inputFile).getChannel();
                long bytesToCopy = input.size();
                if (curFileNum++ < lastFileNumToCopy) {
                    // Compact-reads files end with a delimiter (8 x 0xff)
                    // followed by a 4 byte int 0 (4 x 0x00). Strip
                    // these on all but the last file.
                    bytesToCopy -= (MessageChunksWriter.DELIMITER_LENGTH + 1+ MessageChunksWriter.SIZE_OF_MESSAGE_LENGTH);
                }

                // Copy the file about 10 megabytes at a time. It would probably
                // be marginally faster to just tell NIO to copy the ENTIRE file
                // in one go, but with very large files Java will freeze until the
                // entire chunck is copied so this makes for a more responsive program
                // should you want to ^C in the middle of the copy. Also, with the single
                // transferTo() you might not see any file size changes in the output file
                // until the entire copy is complete.
                long position = 0;
                while (position < bytesToCopy) {
                    long bytesToCopyThisTime = Math.min(maxChunkSize, bytesToCopy - position);
                    position += input.transferTo(position, bytesToCopyThisTime, output);
                }
                input.close();
                input = null;
            }
            System.out.printf("Concatenated %d files.%n", lastFileNumToCopy + 1);
        } finally {
            if (input != null) {
                input.close();
            }
            if (output != null) {
                output.close();
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

    public boolean isQuickConcat() {
        return quickConcat;
    }

    public void setQuickConcat(boolean quickConcat) {
        this.quickConcat = quickConcat;
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
