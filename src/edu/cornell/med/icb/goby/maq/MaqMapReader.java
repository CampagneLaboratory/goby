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

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.goby.maq.MaqStream.LEArrayReader;
import edu.cornell.med.icb.goby.maq.MaqStream.MaqLEDataInputStream;
import edu.cornell.med.icb.goby.maq.filters.AbstractMaqMapEntryFilter;

import java.io.Closeable;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Properties;
import java.util.zip.GZIPInputStream;

/**
 * A reader for MAQ MAP files.
 *
 * Note: By default, this will --NOT-- create a new object for every read but will reuse
 * the same MaqMapEntry object over and over to save time and memory (reduced object creation).
 * If you need to store the MaqMapEntry values that come from this reader in a list, etc.
 * you must set mutableReads to FALSE.
 *
 * @author Kevin Dorff
 */
public class MaqMapReader implements Iterator<MaqMapEntry>, Iterable<MaqMapEntry>, Closeable {

    /**
     * The little endian input stream reader for reading the data.
     */
    protected DataInput is;

    /**
     * The Maq MAP file header (master).
     */
    private MaqMapHeader maqMapHeader;

    /**
     * The Maq MAP file header for the CURRENT FILE.
     */
    private MaqMapHeader currentFileMaqMapHeader;

    /**
     * The next MaqMapEntry entry.
     */
    private MaqMapEntry nextEntry;

    /**
     * Buffer to store the fixed length entry data. We'll read the entire record into this
     * buffer and then use a LEDataInputStream-ByteArrayInputStream to read the data.
     */
    private final byte[] entryData;

    /**
     * The maximum read length for the file we are reading.
     */
    protected final int maxReadLen;

    /**
     * The maximum read length for the file we are reading.
     */
    protected MaqMapEntry mutableEntry;

    /**
     * Use one MaqMapEntry for all reads (don't create a new one for each read), default is TRUE.
     */
    protected boolean mutableReads;

    /** The read name identifiers. */
    protected IndexedIdentifier readNameIdentifiers;

    /** The entry filter to use or null. */
    private AbstractMaqMapEntryFilter entryFilter;

    /** Assists with reading maq records faster. */
    private final LEArrayReader maqRecordReader;

    /**
     * Create a MaqMapReader for reading a Maq MAP file from a filename.
     *
     * @param filenameVal   the filename to read from
     * @param maxReadLenVal the max read length value (64 or 128, generally 128).
     * @param header        the previous MaqMapHeader to append to or null to start
     *                      with a new header
     * @throws IOException error reading
     */
    public MaqMapReader(
            final String filenameVal, final int maxReadLenVal, final MaqMapHeader header)
            throws IOException {
        this(maxReadLenVal, new MaqLEDataInputStream(new GZIPInputStream(
                new FileInputStream(filenameVal))), header);
    }

    /**
     * Create a MaqMapReader for reading a Maq MAP file from a filename.
     *
     * @param is the input stream to read from. This should be a raw input stream
     * as MaqLEDataInputStream and optionally GZIPInputStream will be layered onto the
     * input stream.
     * @param maxReadLenVal the max read length value (64 or 128, generally 128).
     * @param header the previous MaqMapHeader to append to or null to start with a new header
     * @param compressed    True if the stream if GZIP compressed, false otherwise.
     * If you want to reads standard MAQ MAP files you must set this to true.
     * @throws IOException error reading
     */
    public MaqMapReader(
            final DataInputStream is, final int maxReadLenVal, final MaqMapHeader header,
            final boolean compressed) throws IOException {
        this(maxReadLenVal, new MaqLEDataInputStream(compressed ? new GZIPInputStream(is) : is),
                header);
    }

    /**
     * Configure.
     * @param props properties to configure with
     */
    public void configure(final Properties props) {
    }

    /**
     * Create a MaqMapReader for reading a Maq MAP file from an input stream.
     *
     * @param maxReadLenVal the max read length value (64 or 128, generally 128).
     * @param dataInput this input stream should be ready to read the file, ie,
     * should already have MaqLEDataInputStream (and GZIPInputStream) as required by
     * MAQ MAP files.
     * @param header the previous MaqMapHeader to append to or null to start with a
     * new header
     * @throws IOException error reading
     */
    public MaqMapReader(
            final int maxReadLenVal, final DataInput dataInput, final MaqMapHeader header)
            throws IOException {
        is = dataInput;
        this.maxReadLen = maxReadLenVal;
        currentFileMaqMapHeader = readHeader();
        if (header == null) {
            maqMapHeader = currentFileMaqMapHeader;
        } else {
            // Add the data in THIS file's header to the previously made header
            maqMapHeader = header;
            maqMapHeader.merge(currentFileMaqMapHeader);
        }
        nextEntry = null;
        entryData = new byte[this.maxReadLen + 8 + 4 + 4 + 4 + MaqConstants.MAX_NAMELEN];
        maqRecordReader = new LEArrayReader(entryData);
        readNameIdentifiers = new IndexedIdentifier();
        mutableEntry = new MaqMapEntry(this.maxReadLen, readNameIdentifiers);
        entryFilter = null;
        mutableReads = true;
    }

    /**
     * Get readToGeneCounts (may be null).
     * @return the readToGeneCounts.
     */
    public AbstractMaqMapEntryFilter getEntryFilter() {
        return entryFilter;
    }

    /**
     * Set the entry filter (may be null).
     * @param filter the filter to use
     */
    public void setEntryFilter(final AbstractMaqMapEntryFilter filter) {
        this.entryFilter = filter;
        if (filter != null) {
            filter.setHeader(maqMapHeader);
        }
    }

    /**
     * Get the read name identifiers (index <-> read names).
     * @return the read name identifiers
     */
    public IndexedIdentifier getReadNameIdentifiers() {
        return readNameIdentifiers;
    }

    /**
     * Set the read name identifiers (index <-> read names).
     * @param readNameIdentifiers the read name identifiers
     */
    public void setReadNameIdentifiers(final IndexedIdentifier readNameIdentifiers) {
        this.readNameIdentifiers = readNameIdentifiers;
        mutableEntry = new MaqMapEntry(this.maxReadLen, readNameIdentifiers);
    }

    /**
     * Get if one MaqMapEntry should be used for all reads (don't create a new one for each read).
     * @return if MaqMapEntry will be reused on reads
     */
    public boolean isMutableReads() {
        return mutableReads;
    }

    /**
     * Set if one MaqMapEntry should be used for all reads (don't create a new one for each read).
     * The default is TRUE for this which means the MaqMapEntry from this reader, by default,
     * cannot be stored in lists, etc. If you need to store MaqMapEntry values, you must set
     * this to FALSE. With the default value of TRUE, far fewer objects will be created.
     * @param mutableReads if MaqMapEntry will be reused on reads (calls to readEntry())
     */
    public void setMutableReads(final boolean mutableReads) {
        this.mutableReads = mutableReads;
    }

    /**
     * Read a MAQ MAP header from the input stream.
     *
     * @return the MaqMapHeader from the data stream
     * @throws IOException error reading
     */
    protected MaqMapHeader readHeader() throws IOException {
        final MaqMapHeader header = new MaqMapHeader();
        final MaqLEDataInputStream leis;
        if (is instanceof MaqLEDataInputStream) {
            leis = (MaqLEDataInputStream) is;
        } else {
            System.out.println("is has class " + is.getClass().getName());
            throw new IOException("input source isn't a MaqLEDataInputStream");
        }
        header.setFormat(leis.readInt());
        if (header.getFormat() != MaqConstants.MAQMAP_FORMAT_NEW) {
            throw new IOException("Wrong MAQMAP alignment format!");
        }
        final int nRef = is.readInt();
        for (int i = 0; i < nRef; i++) {
            header.addRefName(leis.readVariableLengthString());
        }
        header.referencesDefined();
        // Number of reads, but is often just 0. Cannot be relied upon.
        header.setNumberOfReads(leis.readLong());
        return header;
    }

    /**
     * Read the next MAQ MAP entry from the input stream. This will throw an EOFException if
     * one attempts to read past the end of the file.
     *
     * @return the next MaqMapEntry from the data stream
     * @throws IOException error reading
     */
    protected MaqMapEntry readEntry() throws IOException {
        is.readFully(entryData);
        maqRecordReader.reset();
        final MaqMapEntry entry;
        if (mutableReads) {
            entry = mutableEntry;
        } else {
            entry = new MaqMapEntry(maxReadLen, readNameIdentifiers);
        }
        final short[] seqVal = new short[maxReadLen];
        for (int i = 0; i < maxReadLen; i++) {
            seqVal[i] = maqRecordReader.readUInt8();
        }
        entry.setSeq(seqVal);
        entry.setSize(maqRecordReader.readUInt8());
        entry.setMapQual(maqRecordReader.readUInt8());
        entry.setInfo1(maqRecordReader.readUInt8());
        entry.setInfo2(maqRecordReader.readUInt8());
        final short[] cVal = new short[2];
        cVal[0] = maqRecordReader.readUInt8();
        cVal[1] = maqRecordReader.readUInt8();
        entry.setC(cVal);
        entry.setFlag(maqRecordReader.readUInt8());
        entry.setAltQual(maqRecordReader.readUInt8());

        // Convert the seqId for the local file to the global value across multiple files
        final long currentFileSeqId = maqRecordReader.readUInt32();
        final String refName = currentFileMaqMapHeader.getRefName((int) currentFileSeqId);
        final int globalSeqId = maqMapHeader.getRefNameIndex(refName);
        entry.setSeqId(globalSeqId);

        entry.setPos(maqRecordReader.readUInt32());
        entry.setDist(maqRecordReader.readInt());
        entry.setReadName(maqRecordReader.readFixedLengthString(MaqConstants.MAX_NAMELEN));
        if (entryFilter != null) {
            entryFilter.inspectEntry(entry);
        }
        return entry;
    }

    /**
     * Get the MaqMapHeader.
     *
     * @return the MaqMapHeader
     */
    public MaqMapHeader getMaqMapHeader() {
        return maqMapHeader;
    }

    /**
     * Check if there are more MaqMapEntry to read.
     *
     * @return true if there are more MaqMapEntry to be retrieved.
     */
    public boolean hasNext() {
        if (nextEntry != null) {
            return true;
        }
        try {
            readNextEntry();
            return true;
        } catch (NoSuchElementException e) {
            return false;
        }
    }

    /**
     * Return the next MaqMapEntry from the file. This will throw NoSuchElementException if you
     * read past the end of the file. To avoid this, first check hasNext().
     *
     * @return the next MaqMapEntry
     */
    public MaqMapEntry next() {
        if (nextEntry == null) {
            readNextEntry();
        }
        final MaqMapEntry toReturn = nextEntry;
        nextEntry = null;
        return toReturn;
    }

    /**
     * Close the reader.
     *
     * @throws IOException error closing
     */
    public void close() throws IOException {
        if (is instanceof Closeable) {
            ((Closeable) is).close();
        }
        is = null;
    }

    /**
     * Does nothing. Unsupported.
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Iterator for MaqMapEntry.
     *
     * @return the MaqMapEntry (this class)
     */
    public Iterator<MaqMapEntry> iterator() {
        return this;
    }

    /**
     * Read the next MaqMapEntry into nextEntry. All actual reads occur in ONLY
     * this method. Will throw NoSuchElementException no more MaqMapEntry's were available.
     */
    private void readNextEntry() {
        try {
            while (true) {
                nextEntry = readEntry();
                if (nextEntry != null) {
                    break;
                }
            }
        } catch (EOFException e) {
            throw new NoSuchElementException();
        } catch (IOException e) {
            throw new NoSuchElementException();
        }
    }
}
