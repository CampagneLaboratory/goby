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
import edu.cornell.med.icb.goby.maq.filters.AbstractMaqMapEntryFilter;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * A reader for multiple MAQ MAP files. All of the data from the first file will
 * be returned, then the next file, etc.
 *
 * Note: By default, this will --NOT-- create a new object for every read but will reuse
 * the same MaqMapEntry object over and over to save time and memory (reduced object creation).
 * If you need to store the MaqMapEntry values that come from this reader in a list, etc.
 * you must set mutableReads to FALSE.
 *
 * @author Kevin Dorff
 */
/**
 * Like MaqMapReader but reads from a list of files. All of the data from the first file will
 * be returned, then the next file, etc.
 * @author Kevin Dorff
 */
public class MultiFileMaqMapReader implements Iterator<MaqMapEntry>,
        Iterable<MaqMapEntry>, Closeable {

    /** The files to read from. */
    private final Iterator<File> filesIterator;

    /** The maxReadLen to use for the MaqFileReader. */
    private final int maxReadLen;

    /** The header to use when reading. */
    private MaqMapHeader header;

    /** The "current" reader. */
    private MaqMapReader reader;

    /**
     * The "current" next, populated on hasNext() and cleared on next().
     */
    private MaqMapEntry next;

    /** If the reader should be using mutableReads. False by default. */
    private boolean mutableReads;

    /** The total number of files we are iterating over. */
    private final int totalNumFiles;

    /** The current file number we are iterating over. */
    private int currentFileNum;

    /**
     * The filter to aid with filtering to keep the base MaqMapEntry values.
     * Requires two passes to work.
     */
    private AbstractMaqMapEntryFilter entryFilter;

    /** The read name identifiers. */
    private IndexedIdentifier readNameIdentifiers;

    /** The open filename. */
    private String openFilename;

    /**
     * Constructor.
     *
     * @param filesVal the files to read from
     * @param maxReadLenVal the max read length value (64 or 128, generally 128).
     * @param headerVal        the previous MaqMapHeader to append to or null to start
     *                      with a new header
     * @throws java.io.IOException error reading
     */
    public MultiFileMaqMapReader(
            final List<File> filesVal, final int maxReadLenVal, final MaqMapHeader headerVal) {
        totalNumFiles = filesVal.size();
        currentFileNum = 0;
        filesIterator = filesVal.iterator();
        maxReadLen = maxReadLenVal;
        header = headerVal;
        readNameIdentifiers = new IndexedIdentifier();
        entryFilter = null;
        mutableReads = true;
    }

    /**
     * Close and clear everything.
     * @throws java.io.IOException error closing
     */
    public void close() throws IOException {
        next = null;
        if (reader != null) {
            reader.close();
            reader = null;
            openFilename = null;
        }
    }

    /**
     * Get the entry filter (may be null).
     * @return the entry filter
     */
    public AbstractMaqMapEntryFilter getEntryFilter() {
        return entryFilter;
    }

    /**
     * Set the entry filter (may be null).
     * @param entryFilter the entry filter
     */
    public void setEntryFilter(final AbstractMaqMapEntryFilter entryFilter) {
        this.entryFilter = entryFilter;
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
    }

    /**
     * Get the MaqMapHeader.
     *
     * @return the MaqMapHeader
     */
    public MaqMapHeader getMaqMapHeader() {
        return header;
    }

    /**
     * Detering if there is a "next" element that can be read.
     * @return true if there is a "next" element.
     */
    public boolean hasNext()  {
        try {
            fetchNext();
            return next != null;
        } catch (IOException e) {
            return false;
        }
    }

    /**
     * Obtain the "next" element.
     * @return the "next" element
     * @throws java.util.NoSuchElementException there was not a "next" element
     */
    public MaqMapEntry next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        final MaqMapEntry toReturn = next;
        next = null;
        return toReturn;
    }

    /**
     * Make the "next" entries for each open file have been read.
     * @throws java.io.IOException error reading
     */
    private void fetchNext() throws IOException {
        if (next != null) {
            return;
        }
        if (reader != null) {
            if (!reader.hasNext()) {
                reader.close();
                reader = null;
                openFilename = null;
            }
        }
        while (reader == null) {
            if (filesIterator.hasNext()) {
                currentFileNum++;
                openFilename = filesIterator.next().toString();
                System.out.printf("Opening file to read: [%d/%d] %s%n",
                        currentFileNum, totalNumFiles, openFilename);
                reader = new MaqMapReader(openFilename, maxReadLen, header);
                reader.setMutableReads(mutableReads);
                reader.setReadNameIdentifiers(readNameIdentifiers);
                reader.setEntryFilter(entryFilter);
                header = reader.getMaqMapHeader();
                if (!reader.hasNext()) {
                    reader.close();
                    reader = null;
                    openFilename = null;
                }
            } else {
                // No files left to open
                return;
            }
        }
        next = reader.next();
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
     * @param mutableReadsVal if MaqMapEntry will be reused on reads (calls to readEntry())
     */
    public void setMutableReads(final boolean mutableReadsVal) {
        this.mutableReads = mutableReadsVal;
    }

    /**
     * The current filter number being processed.
     * @return current filter number being processed
     */
    public int getCurrentFileNum() {
        return currentFileNum;
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
     * Does nothing. Unsupported.
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }
}
