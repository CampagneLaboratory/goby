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

package edu.cornell.med.icb.goby.readers;

import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import it.unimi.dsi.io.FastBufferedReader;
import org.apache.commons.io.IOUtils;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

/**
 * A reader for <a href="http://en.wikipedia.org/wiki/FASTA_format">FASTA</a>
 * or <a href="http://en.wikipedia.org/wiki/FASTQ_format">FASTQ</a> files.
 * This reuses the same {@link edu.cornell.med.icb.goby.readers.FastXEntry} over and
 * over for reading the file, so don't directly store it. If you need to store the
 * resultant FastXEntry object, use the {@link edu.cornell.med.icb.goby.readers.FastXEntry#clone()}
 * method to duplicate the FastXEntry object.
 *
 * For FASTQ this parser assumes the # of quality symbols is the same (or more,
 * but should be the same) than the # of sequence characters. The position of linefeeds
 * makes NO DIFFERENCE to this parser in the sequence or quality characters.
 *
 * @author Kevin Dorff
 */
public class FastXReader implements Iterator<FastXEntry>, Iterable<FastXEntry>, Closeable {
    /** The start of line character the designates a new FASTA record. */
    private static final char FASTA_RECORD_START = '>';

    /** The start of line character the designates a new FASTA record. */
    private char fastXRecordStart = FASTA_RECORD_START;

    /** The reader for the FASTA file. */
    private BufferedReader reader;

    /** The "nextEntry" used by next/hasNext. */
    private FastXEntry nextEntry;

    /** The MutableString entry to continually reuse for reading the FASTA file. */
    private final FastXEntry mutableEntry;

    /** The detected file type ("fa" or "fq"). */
    private final String fileType;

    /**
     * Create the FASTX reader.
     * @param file the file that contains the FASTA / FASTQ data
     * @throws IOException error reading or the input stream doesn't support "mark"
     */
    public FastXReader(final String file) throws IOException {
        this(file.endsWith(".gz")
                ? new GZIPInputStream(new FastBufferedInputStream(new FileInputStream(file)))
                : new FileInputStream(file));
    }

    /**
     * Create the FASTX reader.
     * @param is the input stream that contains FASTA / FASTQ data.
     * @throws IOException error reading or the input stream doesn't support "mark"
     */
    public FastXReader(final InputStream is) throws IOException {
        reader = new BufferedReader(new FastBufferedReader(new InputStreamReader(new FastBufferedInputStream(is))));
        nextEntry = null;
        mutableEntry = new FastXEntry();
        if (!reader.markSupported()) {
            reader.close();
            reader = null;
            throw new IOException("FastaReader requires markSupported() on its input stream");
        }
        // We assume the first record will be within the first 32k of the file.
        reader.mark(32768);
        while (true) {
            final String nextLine = reader.readLine();
            if (nextLine == null) {
                // Odd, couldn't find a record. oh well.
                // Assume it is a FASTA file (the default assumption).
                break;
            }
            if (nextLine.length() == 0 || nextLine.charAt(0) == ';' || nextLine.charAt(0) == '#') {
                // Comment or blank line at TOP of file. After this all lines are assumed
                // to be data.
                continue;
            }
            fastXRecordStart = nextLine.charAt(0);
            break;
        }
        reader.reset();

        if (fastXRecordStart == '>') {
            fileType = "fa";
        } else if (fastXRecordStart == '@') {
            fileType = "fq";
        } else {
            fileType = "UNKNOWN";
        }
    }

    /**
     * Get the file type (extension for file). This is determined from the
     * content of the file, not from the file name when the file was opened.
     * @return the file type
     */
    public String getFileType() {
        return fileType;
    }

    /**
     * Check if there are more FASTA entries to read.
     *
     * @return true if there are more FASTA entries to be retrieved.
     */
    public boolean hasNext() {
        if (nextEntry != null) {
            return true;
        }
        try {
            readNextEntry();
            return nextEntry != null;
        } catch (IOException e) {
            return false;
        }
    }

    /**
     * Return the next FASTA entries. This returns next FastaEntry that is
     * re-used for every next! Do not directly store this value as it is reused.
     * @return the next FastaEntry entry
     */
    public FastXEntry next() {
        if (nextEntry == null) {
            try {
                readNextEntry();
                if (nextEntry == null) {
                    throw new NoSuchElementException();
                }
            } catch (IOException e) {
                throw new GobyRuntimeException(e);
            }
        }
        final FastXEntry toReturn = nextEntry;
        nextEntry = null;
        return toReturn;
    }

    /**
     * Read the next FASTA record.
     * @throws IOException error reading
     */
    private void readNextEntry() throws IOException {
        nextEntry = mutableEntry;
        nextEntry.reset();
        while (true) {
            reader.mark(32768);  // 80 chars is recommended for FASTA, but...
            final String nextLine = reader.readLine();
            if (nextLine == null) {
                // No more lines to read
                if (nextEntry.getEntry().length() == 0) {
                    nextEntry = null;
                } else {
                    // End of file, it IS complete.
                    nextEntry.setEntryComplete(true);
                }
                return;
            }
            if (nextLine.length() == 0) {
                // blank line
                continue;
            }
            final char firstChar = nextLine.charAt(0);
            if ((firstChar == ';' || firstChar == '#') && (nextEntry.getEntry().length() == 0)) {
                // We are between entries and found a comment
                continue;
            }
            if (!nextEntry.addLine(nextLine)) {
                reader.reset();
                return;
            }
            if (nextEntry.isEntryComplete()) {
                return;
            }
        }
    }

    /**
     * Does nothing. Unsupported.
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Iterator for FASTA entries.
     *
     * @return the MutableString iterator (this class)
     */
    public Iterator<FastXEntry> iterator() {
        return this;
    }

    /**
     * Close the reader.
     *
     * @throws IOException error closing
     */
    public void close() throws IOException {
        IOUtils.closeQuietly(reader);
    }
}
