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
import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import edu.cornell.med.icb.goby.readers.FastXEntry;
import edu.cornell.med.icb.goby.readers.FastXReader;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.tissueinfo.similarity.GeneTranscriptRelationships;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntCollection;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

/**
 * Class to split transcripts FASTA into multiple files. This class will distribute the transcripts for the same gene
 * across multiple files such that two transcripts from the same gene are never present in the same file. The split
 * files can be used as reference to implement transcript searches.
 *
 * @author Kevin Dorff
 */
public class SplitTranscriptsMode extends AbstractGobyMode {
    /**
     * Used to log debugging and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(SplitTranscriptsMode.class);
    /**
     * The mode name.
     */
    public static final String MODE_NAME = "split-transcripts";

    /**
     * The mode description help text.
     */
    public static final String MODE_DESCRIPTION = "Split transcripts FASTA into multiple files.";

    /**
     * The configuration.
     */
    private SplitTranscriptsConfig config;

    /**
     * Transcript and gene id from header.
     */
    private Map<String, MutableString> transcriptHeader;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Configure via command line.
     *
     * @param args command line args
     * @return this object for chaining
     * @throws java.io.IOException error jsap configuring / jsap parsing
     * @throws com.martiansoftware.jsap.JSAPException error jsap configuring / jsap parsing
     */
    @Override
    public SplitTranscriptsMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        return configure(new SplitTranscriptsConfig(jsapResult));
    }

    /**
     * Configure via CountReadsModeConfig object.
     *
     * @param configVal the config object
     * @return this object for chaining
     */
    public SplitTranscriptsMode configure(final SplitTranscriptsConfig configVal) {
        this.config = configVal;
        transcriptHeader = new HashMap<String, MutableString>();
        return this;
    }

    /**
     * Perform the split transcripts mode.
     *
     * @throws IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        // Load the gene to transcripts file
        if (!config.validate()) {
            throw new IOException("Invalid SplitTranscripts configuration");
        }
        final GeneTranscriptRelationships gtr = new GeneTranscriptRelationships();
        final IndexedIdentifier transcriptIdents = new IndexedIdentifier();
        final Int2ObjectMap<MutableString> transcriptIndexToIdMap
                = new Int2ObjectOpenHashMap<MutableString>();
        final List<FastXEntry> fastxEntries = new LinkedList<FastXEntry>();
        //
        // Pass through the file once to collect the transcript - gene relationships
        //
        int entryCount = 0;
        try {
            for (final FastXEntry entry : new FastXReader(config.getInputFile())) {
                entryCount++;
                parseHeader(entry.getEntryHeader());
                final MutableString transcriptId = transcriptHeader.get("transcriptId");
                final MutableString geneId = transcriptHeader.get("geneId");

                final int transcriptIndex = transcriptIdents.registerIdentifier(transcriptId);
                gtr.addRelationship(geneId, transcriptIndex);

                transcriptIndexToIdMap.put(transcriptIndex, transcriptId);

                fastxEntries.add(entry.clone());
            }
        } catch (CloneNotSupportedException e) {
            LOG.error("Couldn't clone for some reason", e);
            throw new GobyRuntimeException("Couldn't clone for some reason", e);
        }

        LOG.info("Loading map of genes-transcripts complete.");

        //
        // Scan through the transcript-gene relationships to determine which
        // transcript id goes into which file
        //
        final Int2IntMap transcriptIndex2FileIndex = new Int2IntOpenHashMap();
        final String configOutputFilename = config.getOutputBase() + ".config";
        final String configOutputPath = FilenameUtils.getFullPath(configOutputFilename);
        if (StringUtils.isNotBlank(configOutputPath)) {
            LOG.info("Creating output directory: " + configOutputPath);
            FileUtils.forceMkdir(new File(configOutputPath));
        }

        PrintWriter configOutput = null;
        try {
            configOutput = new PrintWriter(configOutputFilename);
            configOutput.println("Ensembl Gene ID\tEnsembl Transcript ID");

            final Int2IntMap fileIndex2NumberOfEntries = new Int2IntOpenHashMap();
            fileIndex2NumberOfEntries.defaultReturnValue(0);
            transcriptIndex2FileIndex.defaultReturnValue(-1);

            final int initialNumberOfFiles = getNumberOfFiles(gtr, transcriptIndex2FileIndex);

            for (int geneIndex = 0; geneIndex < gtr.getNumberOfGenes(); geneIndex++) {
                final MutableString geneId = gtr.getGeneId(geneIndex);
                final IntSet transcriptIndices = gtr.getTranscriptSet(geneIndex);
                int fileNum = 0;

                for (final int transcriptIndex : transcriptIndices) {
                    if (transcriptIndex2FileIndex.get(transcriptIndex) != -1) {
                        LOG.warn("Skipping repeated transcriptIndex: " + transcriptIndex);
                        continue;
                    }
                    final int maxEntriesPerFile = config.getMaxEntriesPerFile();
                    final int numberOfEntriesInOriginalBucket = fileIndex2NumberOfEntries.get(fileNum);
                    final int adjustedFileIndex = fileNum + initialNumberOfFiles * (numberOfEntriesInOriginalBucket / maxEntriesPerFile);

                    transcriptIndex2FileIndex.put(transcriptIndex, adjustedFileIndex);
                    fileIndex2NumberOfEntries.put(fileNum, fileIndex2NumberOfEntries.get(fileNum) + 1);
                    final MutableString transcriptId = transcriptIndexToIdMap.get(transcriptIndex);
                    configOutput.printf("%s\t%s%n", geneId, transcriptId);

                    fileNum++;
                }
            }
        } finally {
            IOUtils.closeQuietly(configOutput);
        }

        final int numFiles = getFileIndices(transcriptIndex2FileIndex).size();
        if (LOG.isInfoEnabled()) {
            LOG.info(NumberFormat.getInstance().format(entryCount)
                    + " entries will be written to " + numFiles + " files");
            final int maxEntriesPerFile = config.getMaxEntriesPerFile();
            if (maxEntriesPerFile < Integer.MAX_VALUE) {
                LOG.info("Each file will contain at most " + maxEntriesPerFile + " entries");
            }
        }

        // formatter for uniquely numbering files each with the same number of digits
        final NumberFormat fileNumberFormatter = getNumberFormatter(numFiles - 1);

        final ProgressLogger progressLogger = new ProgressLogger();
        progressLogger.expectedUpdates = entryCount;
        progressLogger.itemsName = "entries";
        progressLogger.start();

        // Write each file one at a time rather than in the order they appear in the input file
        // to avoid the issue of having too many streams open at the same or continually opening
        // and closing streams which is quite costly.  We could store the gene/transcripts in
        // memory and then just write the files at the end but that could be worse.
        for (final int fileIndex : getFileIndices(transcriptIndex2FileIndex)) {
            final String filename =
                    config.getOutputBase() + "." + fileNumberFormatter.format(fileIndex) + ".fa.gz";
            PrintStream printStream = null;
            try {
                // each file is compressed
                printStream = new PrintStream(new GZIPOutputStream(new FileOutputStream(filename)));

                //
                // Read through the input file get the actual sequence information
                //
                final Iterator<FastXEntry> entries = fastxEntries.iterator();
                while (entries.hasNext()) {
                    final FastXEntry entry = entries.next();
                    parseHeader(entry.getEntryHeader());
                    final MutableString transcriptId = transcriptHeader.get("transcriptId");
                    final MutableString geneId = transcriptHeader.get("geneId");
                    final int transcriptIndex = transcriptIdents.getInt(transcriptId);
                    final int transcriptFileIndex = transcriptIndex2FileIndex.get(transcriptIndex);
                    if (transcriptFileIndex == fileIndex) {
                        printStream.print(entry.getHeaderSymbol());
                        printStream.print(transcriptId);
                        printStream.print(" gene:");
                        printStream.println(geneId);
                        printStream.println(entry.getEntrySansHeader());
                        entries.remove();
                        progressLogger.lightUpdate();
                    }
                }
            } finally {
                IOUtils.closeQuietly(printStream);
            }
        }

        assert progressLogger.count == entryCount : "Some entries were not processed!";
        progressLogger.done();
    }

    private IntCollection getFileIndices(final Int2IntMap transcriptIndex2FileIndex) {
        final IntSet result = new IntArraySet();
        result.addAll(transcriptIndex2FileIndex.values());
        return result;
    }

    /**
     * Calcuate the number of files the gene to transcript relationships should be split to.
     * @param gtr The gene to transcript relationships
     * @param transcriptIndex2FileIndex A map of transcript identifiers to file indices
     * @return The number of files to write
     */
    private int getNumberOfFiles(final GeneTranscriptRelationships gtr,
                                 final Int2IntMap transcriptIndex2FileIndex) {
        int numFiles = 0;
        for (int geneIndex = 0; geneIndex < gtr.getNumberOfGenes(); geneIndex++) {
            final IntSet transcriptIndices = gtr.getTranscriptSet(geneIndex);
            int fileNum = 0;
            for (final int transcriptIndex : transcriptIndices) {
                if (transcriptIndex2FileIndex.get(transcriptIndex) != -1) {
                    LOG.warn("Skipping repeated transcriptIndex: " + transcriptIndex);
                    continue;
                }

                numFiles = Math.max(fileNum, numFiles);
                fileNum++;
            }
        }
        return ++numFiles;
    }

    /**
     * Extract the transcript and gene ids from the given string.
     * @param header The string to extract the information from.  Generally speaking
     * this is the comment line from a FASTA entry (without the ">" character)
     */
    private void parseHeader(final MutableString header) {
        final int endOfTranscriptId = header.indexOf(' ');
        transcriptHeader.put("transcriptId", header.substring(0, endOfTranscriptId));

        final int startOfGeneId = header.lastIndexOf(' ');
        transcriptHeader.put("geneId", header.substring(startOfGeneId + 6));
    }

    /**
     * Get a number formatter to print leading zeros up to n.
     *
     * @param n The largest number that will be formatted
     * @return the NumberFormat for n
     */
    public NumberFormat getNumberFormatter(final int n) {
        assert n >= 0 : "n must be non-negative";
        final int numDigits;
        if (n == 0) {
            numDigits = 1;
        } else {
            numDigits = 1 + (int) (Math.log10(n));
        }

        final NumberFormat numberFormat = NumberFormat.getInstance();
        numberFormat.setGroupingUsed(false);
        numberFormat.setMinimumIntegerDigits(numDigits);
        return numberFormat;
    }


    /**
     * Main method.
     *
     * @param args command line args.
     * @throws JSAPException error parsing
     * @throws IOException error parsing or executing.
     */
    public static void main(final String[] args) throws IOException, JSAPException {
        new SplitTranscriptsMode().configure(args).execute();
    }

    /**
     * Config class for SplitTranscriptsMode.
     *
     * @author Kevin Dorff
     */
    private static class SplitTranscriptsConfig {
        /** The fasta input file. */
        private final String inputFile;
        /** The outputBase. */
        private final String outputBase;
        /** Maximum number of entries that should be written to any given output file. */
        private final int maxEntriesPerFile;

        /**
         * Create a config object based on parsed JSAP.
         * @param jsapResult the parsed JSAP
         */
        public SplitTranscriptsConfig(final JSAPResult jsapResult)  {
            inputFile = jsapResult.getString("input");
            outputBase = jsapResult.getString("output");
            maxEntriesPerFile = jsapResult.getInt("max-entries-per-file");
        }

        /**
         * Get the input file to process.
         * @return the input filename
         */
        public String getInputFile() {
            return inputFile;
        }

        /**
         * Get the basename/directory for the output file(s).
         * @return the basename
         */
        public String getOutputBase() {
            return outputBase;
        }

        /**
         * Get the maximum number of entries that should be written to any given output file.
         * @return the maximum number of entries to write
         */
        public int getMaxEntriesPerFile() {
            return maxEntriesPerFile;
        }

        /**
         * Validate the values.
         * @return true of values validated
         */
        public boolean validate() {
            boolean validates = true;
            if (!(new File(inputFile)).exists()) {
                LOG.error("input file " + inputFile + " does not exist.");
                validates = false;
            }
            return validates;
        }
    }
}
