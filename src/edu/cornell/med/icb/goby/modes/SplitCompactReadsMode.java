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
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.lang.StringUtils;

import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Converts a Compact file to Fasta format.
 *
 * @author Fabien Campagne
 *         Date: May 4 2009
 *         Time: 12:28 PM
 */
public class SplitCompactReadsMode extends AbstractGobyMode {
    /**
     * An "unset" value for startPosition and endPosition.
     */
    private static final long UNSET_POSITION = -1;

    /** The input filename. */
    private String inputFilename;

    /** The output filename. */
    private String outputFilename;

    /** The start position. */
    private long startPosition = UNSET_POSITION;
    /** The end position. */
    private long endPosition = UNSET_POSITION;

    /** sequences per chunk in the written file. */
    private int sequencePerChunk = 10000;

    /** The mode name. */
    public static final String MODE_NAME = "split-compact-reads";
    public static final String MODE_DESCRIPTION = "Converts a Compact file to Fasta format.";

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
     * @throws IOException error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilename = jsapResult.getString("input");
        outputFilename = jsapResult.getString("output");
        startPosition = jsapResult.getLong("start-position");
        endPosition = jsapResult.getLong("end-position");
        sequencePerChunk = jsapResult.getInt("sequence-per-chunk");
        return this;
    }

    /**
     * Actually perform the split of the compact reads file between
     * start and end position to the new compact reads file.
     * @throws IOException
     */
    @Override
    public void execute() throws IOException {
        if (StringUtils.isBlank(inputFilename)) {
            throw new IOException("--input not specified");
        }
        if (StringUtils.isBlank(outputFilename)) {
            throw new IOException("--output not specified");
        }

        if (startPosition == UNSET_POSITION || startPosition < 0) {
            throw new IOException("--start-position not specified");
        }
        if (endPosition == UNSET_POSITION || endPosition < 0) {
            throw new IOException("--end-position not specified");
        }
        if (endPosition <= startPosition) {
            throw new IOException(
                    "--end-position cannot be less than or the same as --start-position");
        }

        final ReadsWriter writer = new ReadsWriter(new FileOutputStream(outputFilename));
        writer.setNumEntriesPerChunk(sequencePerChunk);
        final MutableString sequence = new MutableString();
        final ReadsReader readsReader = new ReadsReader(startPosition, endPosition, inputFilename);

        try {
            for (final Reads.ReadEntry readEntry : readsReader) {
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
                if (readEntry.hasSequence()) {
                    ReadsReader.decodeSequence(readEntry, sequence);
                    writer.setSequence(sequence);
                } else {
                    writer.setSequence("");
                }
                writer.appendEntry();
            }
            writer.printStats(System.out);
        } finally {
            writer.close();
            readsReader.close();
        }
    }

    /**
     * Set the end position. This will stop copying records
     * ending at the endPosition. If endPosition is in the
     * middle of a record, this will copy to the end of that record.
     * @param endPosition the start position
     */
    public void setEndPosition(final long endPosition) {
        this.endPosition = endPosition;
    }

    /**
     * Get the end position. This will stop copying records
     * ending at the endPosition. If endPosition is in the
     * middle of a record, this will copy to the end of that record.
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
     * @param startPosition the start position
     */
    public void setStartPosition(final long startPosition) {
        this.startPosition = startPosition;
    }

    /**
     * Get the start position. This will start copying records
     * starting at the startPosition. If startPosition is in the
     * middle of a record, this will advance to the start of the
     * next record.
     * @return the start position
     */
    public long getStartPosition() {
        return startPosition;
    }

    /**
     * Set the output filename.
     * @param inputFilename the output filename
     */
    public void setInputFilename(final String inputFilename) {
        this.inputFilename = inputFilename;
    }

    /**
     * Get the input filename.
     * @return the input filename
     */
    public String getInputFilename() {
        return inputFilename;
    }

    /**
     * Set the output filename.
     * @param outputFilename the output filename
     */
    public void setOutputFilename(final String outputFilename) {
        this.outputFilename = outputFilename;
    }

    /**
     * Get the output filename.
     * @return the output filename
     */
    public String getOutputFilename() {
        return outputFilename;
    }

    /**
     * Main mode for splitting compact reads files from a start position
     * to and end position.
     * @param args command line arguments
     * @throws IOException IO error
     * @throws JSAPException command line parsing error.
     */
    public static void main(final String[] args) throws IOException, JSAPException {
        new SplitCompactReadsMode().configure(args).execute();
    }
}
