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

import cern.jet.random.engine.MersenneTwister;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import edu.cornell.med.icb.util.RandomAdapter;
import it.unimi.dsi.fastutil.chars.CharArraySet;
import it.unimi.dsi.fastutil.chars.CharSet;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Date;

/**
 * Reformat a compact file, possibly dropping identifiers, or descriptions, or splitting the file.
 * When a compact-reads file is split, reads in each split are renumbered (their readIndex
 * is changed), starting at zero for the first sequence of each split. This ensures that indices
 * are correctly concatenated back together.
 *
 * @author Fabien Campagne
 *         Date: Apr 28, 2009
 *         Time: 6:03:56 PM
 */
public class ReformatCompactReadsMode extends AbstractGobyMode {
    private String[] inputFilenames;

    private boolean pushDescription;
    private boolean pushIdentifier;
    private String outputFile;
    private boolean mutateSequences;
    private int numberOfMismatches;
    private CharSet bases;

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "reformat-compact-reads";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Reformat a compact file, possibly dropping "
            + "identifiers, or descriptions, or splitting the file. The mode also supports randomly mutating "
            + "bases and other modifications, trimming reads, etc. When a compact-reads file "
            + "is split, reads in each split are renumbered (their readIndex is changed), "
            + "starting at zero for the first sequence of each split. This ensures that indices "
            + "are correctly concatenated back together. ";

    private int sequencePerChunk = 10000;
    private boolean excludeSequences;
    private int sequencePerOutput = Integer.MAX_VALUE;
    /**
     * Any sequences shorter than this length will be excluded from the output.
     */
    private int minReadLength;
    /**
     * Any sequences longer than this length will be excluded from the output.
     */
    private int maxReadLength = Integer.MAX_VALUE;
    /**
     * Any sequences longer than this length will be cut to this length in the output.
     */
    private int trimReadLength = Integer.MAX_VALUE;
    private final ObjectSet<String> outputFilenames = new ObjectOpenHashSet<String>();

    /**
     * An "unset" value for startPosition and endPosition.
     */
    private boolean hasStartOrEndPosition;

    /**
     * The start position for the reformat.
     */
    private long startPosition;

    /**
     * The end position for the reformat.
     */
    private long endPosition = Long.MAX_VALUE;
    /**
     * The number of bases to trim at the start of the sequence.
     */
    private int trimReadStartLength;

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
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        inputFilenames = jsapResult.getStringArray("input");
        pushDescription = jsapResult.getBoolean("include-descriptions");
        pushIdentifier = jsapResult.getBoolean("include-identifiers");
        excludeSequences = jsapResult.getBoolean("exclude-sequences");
        outputFile = jsapResult.getString("output");
        sequencePerChunk = jsapResult.getInt("sequence-per-chunk");
        sequencePerOutput = jsapResult.getInt("sequence-per-output", sequencePerOutput);
        minReadLength = jsapResult.getInt("minimum-read-length", minReadLength);
        maxReadLength = jsapResult.getInt("maximum-read-length", maxReadLength);
        trimReadLength = jsapResult.getInt("trim-read-length", trimReadLength);
        trimReadStartLength = jsapResult.getInt("trim-read-start", trimReadLength);
        mutateSequences = jsapResult.getBoolean("mutate-sequences");
        numberOfMismatches = jsapResult.getInt("mismatch-number");

        if (jsapResult.contains("start-position") || jsapResult.contains("end-position")) {
            hasStartOrEndPosition = true;
            startPosition = jsapResult.getLong("start-position", 0L);
            endPosition = jsapResult.getLong("end-position", Long.MAX_VALUE);
        }

        if (startPosition < 0L) {
            throw new JSAPException("Start position must not be less than zero");
        }
        if (endPosition < 0L) {
            throw new JSAPException("End position must not be less than zero");
        }
        if (startPosition > endPosition) {
            throw new JSAPException("Start position must not be greater than the end position");
        }

        bases = new CharArraySet();
        return this;
    }

    /**
     * Reformat compact reads.
     *
     * @throws IOException
     */
    @Override
    public void execute() throws IOException {
        final int numToProcess = inputFilenames.length;
        int numProcessed = 0;
        final MutableString sequence = new MutableString();
        for (final String inputFilename : inputFilenames) {
            int splitIndex = 0;
            final String outputBasename;
            final String outputFilename;

            // if there is only one input file to process
            if (numToProcess == 1 && StringUtils.isNotBlank(outputFile)) {
                outputBasename = outputFile;

                boolean hasCompactExtension = false;
                for (final String extension : FileExtensionHelper.COMPACT_READS_FILE_EXTS) {
                    if (outputFile.endsWith(extension)) {
                        hasCompactExtension = true;
                        break;
                    }
                }

                // if it already has a compact-reads extension leave it alone
                // otherwise build it based on the basename
                if (hasCompactExtension) {
                    outputFilename = outputFile;
                } else {
                    outputFilename = getOutputFilename(outputBasename, ++splitIndex);
                }
            } else {
                // there are a bunch of files to process so number them uniquely based on the name
                outputBasename = stripCompactReadExtensions(inputFilename);
                outputFilename = getOutputFilename(outputBasename, ++splitIndex);
            }

            if (inputFilename.equals(outputFilename)) {
                System.err.println("input cannot equal the output name");
                System.exit(2);
            }
            System.out.printf("Converting [%d/%d] %s to %s%n",
                    ++numProcessed, numToProcess, inputFilename, outputFilename);

            outputFilenames.add(outputFilename);
            ReadsWriter writer = new ReadsWriter(new FileOutputStream(outputFilename));
            writer.setNumEntriesPerChunk(sequencePerChunk);

            final ReadsReader readsReader;
            final FastBufferedInputStream inputFileStream =
                    new FastBufferedInputStream(new FileInputStream(inputFilename));

            if (hasStartOrEndPosition) {
                readsReader = new ReadsReader(startPosition, endPosition, inputFileStream);
            } else {
                readsReader = new ReadsReader(inputFileStream);
            }

            int entriesInOutputFile = 0;
            for (final Reads.ReadEntry entry : readsReader) {
                final int readLength = entry.getReadLength();
                if (readLength < minReadLength || readLength > maxReadLength) {
                    continue;
                }

                if (pushDescription && entry.hasDescription()) {
                    writer.setDescription(entry.getDescription());
                }

                if (pushIdentifier && entry.hasReadIdentifier()) {
                    writer.setIdentifier(entry.getReadIdentifier());
                }

                if (!excludeSequences) {
                    ReadsReader.decodeSequence(entry, sequence);
                    // trim the read length down to size
                    if (sequence.length() > trimReadLength) {
                        sequence.length(trimReadLength);
                    }
                    if (trimReadStartLength > 0) {
                        sequence.delete(0, trimReadStartLength);
                    }
                    if (mutateSequences) {
                        mutate(sequence, numberOfMismatches);
                    }
                    writer.setSequence(sequence);
                } else {
                    writer.setSequence("");

                }

                if (entry.hasQualityScores()) {
                    writer.setQualityScores(entry.getQualityScores().toByteArray());
                }
                // Important: preserve the read index in the input entry:
                writer.appendEntry(entry.getReadIndex());
                //writer.appendEntry();
                entriesInOutputFile++;

                if (entriesInOutputFile > sequencePerOutput) {
                    writer.close();

                    final String newOutputFilename = getOutputFilename(outputBasename, ++splitIndex);
                    outputFilenames.add(newOutputFilename);
                    System.out.printf("Splitting into %s%n", newOutputFilename);
                    writer = new ReadsWriter(new FileOutputStream(newOutputFilename));
                    writer.setNumEntriesPerChunk(sequencePerChunk);
                    entriesInOutputFile = 0;
                }
            }

            writer.close();
            writer.printStats(System.out);
        }
    }

    /**
     * Introduce the given number of point mutations in the given sequence.
     *
     * @param sequence
     * @param numberOfMismatches
     */
    private void mutate(final MutableString sequence, final int numberOfMismatches) {
        for (final char base : sequence.toCharArray()) {
            bases.add(base);
        }
        if (sequence.length() < numberOfMismatches) {
            System.err.printf("Cannot introduce %d mismatches in a sequence that has only %d "
                    + "residues. Skipping this sequence.", numberOfMismatches, sequence.length());
            return;
        }
        final IntSet alreadyMutated = new IntArraySet();
        final RandomAdapter random = new RandomAdapter(new MersenneTwister(new Date()));
        for (int i = 0; i < numberOfMismatches; i++) {
            int mutationPosition;
            do {
                mutationPosition = random.choose(0, sequence.length() - 1);
            }
            while (alreadyMutated.contains(mutationPosition));
            char newBase;
            final char oldBase = sequence.charAt(mutationPosition);
            do {
                newBase = bases.toCharArray()[random.choose(0, bases.size() - 1)];
                // introduce the mutation:
                sequence.charAt(mutationPosition, newBase);
            } while (newBase == oldBase);
            alreadyMutated.add(mutationPosition);
        }
    }

    private String getOutputFilename(final String outputBasename, final int splitIndex) {
        final boolean neverSplits = sequencePerOutput == Integer.MAX_VALUE;
        return outputBasename + (neverSplits ? "" : "-" + splitIndex) + ".compact-reads";
    }

    /**
     * Get the filename including path WITHOUT fastx extensions (including .gz if it is there).
     *
     * @param name the full path to the file in question
     * @return the filename without the fastx/gz extensions or the same name of those extensions
     *         weren't found.
     */
    private static String stripCompactReadExtensions(final String name) {
        final String filename = FilenameUtils.getName(name);
        for (final String ext : FileExtensionHelper.COMPACT_READS_FILE_EXTS) {
            if (filename.endsWith(ext)) {
                return FilenameUtils.getFullPath(name)
                        + filename.substring(0, filename.lastIndexOf(ext));
            }
        }
        return name;
    }

    public void setInputFilenames(final String... filenames) {
        inputFilenames = filenames;
    }

    public void setInputFilenames(final File... files) {
        inputFilenames = new String[files.length];
        int i = 0;
        for (final File in : files) {
            inputFilenames[i++] = in.getPath();
        }
    }

    /**
     * Set the number of sequences per output file.
     *
     * @param sequencePerOutput number of sequences per output file.
     */
    public void setSequencePerOutput(final int sequencePerOutput) {
        this.sequencePerOutput = sequencePerOutput;
    }

    /**
     * Returns the set of filenames where data was output.
     *
     * @return The list of filenames that were created as a result of the reformat operation
     */
    public String[] getOutputFilenames() {
        return outputFilenames.toArray(new String[outputFilenames.size()]);
    }

    public void setOutputFile(final String output) {
        outputFile = output;
    }

    /**
     * Set the end position. This will stop copying records
     * ending at the endPosition. If endPosition is in the
     * middle of a record, this will copy to the end of that record.
     *
     * @param endPosition the start position
     */
    public void setEndPosition(final long endPosition) {
        assert endPosition >= 0 : "End position must not be negative";
        this.endPosition = endPosition;
        hasStartOrEndPosition = true;
    }

    /**
     * Get the end position. This will stop copying records
     * ending at the endPosition. If endPosition is in the
     * middle of a record, this will copy to the end of that record.
     *
     * @return the start position
     */
    public long getEndPosition() {
        return endPosition;
    }

    /**
     * Set the start position. This will start copying records
     * starting at the startPosition. If startPosition is in the
     * middle of a record, this will advance to the start of the
     * next record.
     *
     * @param startPosition the start position
     */
    public void setStartPosition(final long startPosition) {
        assert startPosition >= 0 : "Start position must not be negative";
        this.startPosition = startPosition;
        hasStartOrEndPosition = true;
    }

    /**
     * Get the start position. This will start copying records
     * starting at the startPosition. If startPosition is in the
     * middle of a record, this will advance to the start of the
     * next record.
     *
     * @return the start position
     */
    public long getStartPosition() {
        return startPosition;
    }

    public int getMinReadLength() {
        return minReadLength;
    }

    public void setMinReadLength(final int minReadLength) {
        this.minReadLength = minReadLength;
    }

    public int getMaxReadLength() {
        return maxReadLength;
    }

    public void setMaxReadLength(final int maxReadLength) {
        this.maxReadLength = maxReadLength;
    }

    public int getTrimReadLength() {
        return trimReadLength;
    }

    public void setTrimReadLength(final int trimReadLength) {
        this.trimReadLength = trimReadLength;
    }

    public int getSequencePerChunk() {
        return sequencePerChunk;
    }

    public void setSequencePerChunk(final int sequencePerChunk) {
        this.sequencePerChunk = sequencePerChunk;
    }

    public static void main(final String[] args) throws IOException, JSAPException {
        new ReformatCompactReadsMode().configure(args).execute();
    }
}
