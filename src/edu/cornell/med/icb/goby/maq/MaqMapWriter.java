/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.goby.maq;

import edu.cornell.med.icb.goby.maq.MaqStream.MaqLEDataOutputStream;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Class to help with writing the little-endian MAQ MAP files. this class actually first
 * writes all the entries to their own temporary file, the writes the header to the actual
 * output file, then append the entries temporary file to the primary output file and finally
 * delete the entries temporary file. This way the header doesn't need to be complete until
 * we are done writing the entries. The numberOfReads written to the header of this file
 * will exactly match the number of entries written.
 *
 * @author Kevin Dorff
 */
public class MaqMapWriter implements Closeable {

    /** Buffer size for writing. */
    private static final int BUFFER_SIZE = 16 * 1024 * 1024;

    /**
     * The Little Endian Data Input Stream we use to read the MAQ data.
     */
    private MaqLEDataOutputStream primaryDos;

    /**
     * The temporary file we use to write entries.
     */
    private File entriesTempFile;

    /**
     * The output stream we use to write entries, goes to a temporary file.
     */
    private MaqLEDataOutputStream entriesDos;

    /**
     * The maximum read length for the file we are reading.
     */
    private final int maxReadLen;

    /**
     * If the header has been written.
     */
    private final MaqMapHeader header;

    /**
     * The number of entries written.
     */
    private long numWritten;

    /**
     * Writer to help with writing MAQ MAP files. The numberOfReads in the file
     * that is written matches the exact number of entries that are written by this
     * writer. header.numberOfReads is not altered by this writer.
     *
     * @param filenameVal   the filename to write the maq map to
     * @param maxReadLenVal the maximum read length for the file
     * @param headerVal     the MaqMapHeader to use for writing
     * @throws java.io.IOException error writing to file (notably from GZIPInputStream).
     */
    public MaqMapWriter(
            final String filenameVal, final int maxReadLenVal, final MaqMapHeader headerVal)
            throws IOException {
        this(maxReadLenVal, new MaqLEDataOutputStream(new GZIPOutputStream(
                new FileOutputStream(filenameVal))), headerVal);
    }

    /**
     * Writer to help with writing MAQ MAP files. The numberOfReads in the file
     * that is written matches the exact number of entries that are written by this
     * writer. header.numberOfReads is not altered by this writer.
     *
     * @param os            the output stream to write the maq map to. This should be a raw output
     *                      stream as MaqLEDataOutputStream and optionally GZIPOutputStream will be layered
     *                      onto this output stream for you to make a MAQ MAP compliant file.
     * @param maxReadLenVal the maximum read length for the file
     * @param headerVal     the MaqMapHeader to use for writing
     * @param compressed    set to true if you want the output stream to be GZip compressed.
     *                      NOTE: If you want the resultant files to be readable with other MAQ tools, you MUST
     *                      set this to true and GZip the output contents.
     * @throws java.io.IOException error writing to file (notably from GZIPInputStream).
     */
    public MaqMapWriter(
            final OutputStream os, final int maxReadLenVal, final MaqMapHeader headerVal,
            final boolean compressed) throws IOException {
        this(maxReadLenVal, new MaqLEDataOutputStream(new BufferedOutputStream(
                compressed ? new GZIPOutputStream(os) : os, BUFFER_SIZE)), headerVal);
    }

    /**
     * Writer to help with writing MAQ MAP files. The numberOfReads in the file
     * that is written matches the exact number of entries that are written by this
     * writer. header.numberOfReads is not altered by this writer.
     *
     * @param maxReadLenVal the maximum read length for the file
     * @param dataOutput    the MaqLEDataOutputStream to write to. This should be a "ready to go"
     *                      output stream, no additional layers (such as GZip) will be added to this stream.
     * @param headerVal     the MaqMapHeader to use for writing
     * @throws java.io.IOException error writing to file (notably from GZIPInputStream).
     */
    public MaqMapWriter(
            final int maxReadLenVal, final MaqLEDataOutputStream dataOutput,
            final MaqMapHeader headerVal) throws IOException {
        this.header = headerVal;
        numWritten = 0;
        primaryDos = dataOutput;
        entriesTempFile = File.createTempFile("maq-map-entries-temp-", ".tmp");
        entriesDos = new MaqLEDataOutputStream(new BufferedOutputStream(
                new FileOutputStream(entriesTempFile), BUFFER_SIZE));
        this.maxReadLen = maxReadLenVal;
    }

    /**
     * Close the writer.
     *
     * @throws IOException error closing
     */
    public void close() throws IOException {
        if (entriesDos != null) {
            entriesDos.close();
            entriesDos = null;
        }
        if (primaryDos != null) {
            writeMaqMapHeader();
            copyEntriesFromTemp();
            primaryDos.close();
            primaryDos = null;
            FileUtils.deleteQuietly(entriesTempFile);
        }
    }

    /**
     * Part of {@link #close()}, this will copy the entries from the entries temp file to the
     * primary output file.
     *
     * @throws IOException error reading (or closing).
     */
    private void copyEntriesFromTemp() throws IOException {
        InputStream is = null;
        try {
            is = new BufferedInputStream(new FileInputStream(entriesTempFile), BUFFER_SIZE);
            final byte[] buffer = new byte[BUFFER_SIZE];
            while (true) {
                final int size = is.read(buffer);
                if (size == -1) {
                    break;
                }
                primaryDos.write(buffer, 0, size);
            }
        } finally {
            IOUtils.closeQuietly(is);
        }
    }

    /**
     * Read a MAQ MAP header from the input stream. This needs to be called before any
     * calls to writeMaqMapEntry().
     *
     * @throws java.io.IOException error writing
     */
    private void writeMaqMapHeader() throws IOException {

        if (header.getFormat() != MaqConstants.MAQMAP_FORMAT_NEW) {
            throw new IOException(String.format(
                    "Wrong alignment format! found %d but expected %d",
                    header.getFormat(), MaqConstants.MAQMAP_FORMAT_NEW));
        }
        primaryDos.writeInt(header.getFormat());
        primaryDos.writeInt(header.getNRef());
        header.referencesDefined();
        for (int i = 0; i < header.getNRef(); i++) {
            primaryDos.writeVariableLengthString(header.getRefName(i));
        }
        // Number of record written
        primaryDos.writeLong(numWritten);
    }

    /**
     * Read the next MAQ MAP entry from the input stream. This will throw an EOFException if
     * one attempts to read past the end of the file.
     *
     * @param entry the MaqMapEntry to write
     * @throws java.io.IOException error writing
     */
    public void writeMaqMapEntry(final MaqMapEntry entry) throws IOException {
        numWritten++;
        final short[] seqVal = entry.getSeq();
        for (int i = 0; i < maxReadLen; i++) {
            entriesDos.writeUInt8(seqVal[i]);
        }

        entriesDos.writeUInt8(entry.getSize());
        entriesDos.writeUInt8(entry.getMapQual());
        entriesDos.writeUInt8(entry.getInfo1());
        entriesDos.writeUInt8(entry.getInfo2());

        final short[] cVal = entry.getC();
        entriesDos.writeUInt8(cVal[0]);
        entriesDos.writeUInt8(cVal[1]);

        entriesDos.writeUInt8(entry.getFlag());
        entriesDos.writeUInt8(entry.getAltQual());
        entriesDos.writeUInt32(entry.getReferenceSequenceId());
        entriesDos.writeUInt32(entry.getPos());
        entriesDos.writeInt(entry.getDist());
        entriesDos.writeFixedLengthString(entry.getName().toString(), MaqConstants.MAX_NAMELEN);
    }
}
