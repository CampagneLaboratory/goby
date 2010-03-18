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

package edu.cornell.med.icb.goby.alignments;

import com.google.protobuf.CodedInputStream;
import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import edu.cornell.med.icb.goby.reads.FastBufferedMessageChunksReader;
import edu.cornell.med.icb.goby.reads.MessageChunksReader;
import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Properties;
import java.util.zip.GZIPInputStream;

/**
 * Reads alignments written with AlignmentWriter.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 6:36:04 PM
 */
public class AlignmentReader extends AbstractAlignmentReader {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(AlignmentReader.class);

    private InputStream headerStream;
    private int numberOfAlignedReads;
    private final MessageChunksReader alignmentEntryReader;
    private Alignments.AlignmentCollection collection;
    private Properties stats;
    private String basename;

    public AlignmentReader(final String basename) throws FileNotFoundException {
        super();
        this.basename = basename;
        final FileInputStream stream = new FileInputStream(basename + ".entries");
        alignmentEntryReader = new MessageChunksReader(stream);
        try {
            headerStream = new GZIPInputStream(new FileInputStream(basename + ".header"));
        } catch (IOException e) {
            // try not compressed for compatibility with 1.4-:
            LOG.trace("falling back to legacy 1.4- uncompressed header.");

            headerStream = new FileInputStream(basename + ".header");
        }
        stats = new Properties();
        final File statsFile = new File(basename + ".stats");
        if (statsFile.exists()) {
            Reader statsFileReader = null;
            try {
                statsFileReader = new FileReader(statsFile);
                stats.load(statsFileReader);
            } catch (IOException e) {
                LOG.warn("cannot load properties for basename: " + basename, e);
            } finally {
                IOUtils.closeQuietly(statsFileReader);
            }
        }
    }

    /**
     * Returns the basename for the alignment being read, or null if the basename is unknown.
     * If a path was part of basename provided to the constructor, it is returned.
     *
     * @return basename for the alignment being read
     */
    public String basename() {
        return basename;
    }

    public AlignmentReader(final InputStream entriesStream) {
        super();
        alignmentEntryReader = new MessageChunksReader(entriesStream);
    }

    /**
     * Initialize the reader to read a segment of the input. Sequences represented by a collection
     * which starts between the input position start and end will be returned upon subsequent
     * calls to {@link #hasNext()}, {@link #next()}.
     *
     * @param start  Start offset in the input file.
     * @param end    End offset in the input file.
     * @param stream Stream over the input file
     * @throws IOException If an error occurs reading the input.
     */
    public AlignmentReader(final long start, final long end, final FastBufferedInputStream stream)
            throws IOException {
        super();
        alignmentEntryReader = new FastBufferedMessageChunksReader(start, end, stream);
    }

    private int numberOfEntries() {
        return collection != null ? collection.getAlignmentEntriesCount() : 0;
    }

    /**
     * Returns true if the input has more entries.
     *
     * @return true if the input has more entries, false otherwise.
     */
    public boolean hasNext() {
        final boolean hasNext = alignmentEntryReader.hasNext(collection, numberOfEntries());
        if (!hasNext) {
            collection = null;
        }
        final GZIPInputStream uncompressStream = alignmentEntryReader.getUncompressStream();
        try {
            if (uncompressStream != null) {
                collection = Alignments.AlignmentCollection.parseFrom(uncompressStream);
                if (collection.getAlignmentEntriesCount() == 0) {
                    return false;
                }
            }
        } catch (IOException e) {
            throw new GobyRuntimeException(e);
        }
        return hasNext;
    }

    /**
     * Returns the next alignment entry from the input stream.
     *
     * @return the alignment read entry from the input stream.
     */
    public Alignments.AlignmentEntry next() {
        if (!alignmentEntryReader.hasNext(collection, numberOfEntries())) {
            throw new NoSuchElementException();
        }
        return collection.getAlignmentEntries(alignmentEntryReader.incrementEntryIndex());
    }

    /**
     * This operation is not supported.
     */
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a reader.");
    }

    /**
     * Read the header of this alignment.
     *
     * @throws IOException If an error occurs.
     */
    @Override
    public void readHeader() throws IOException {
        // accept very large header messages, since these may contain query identifiers:
        final CodedInputStream codedInput = CodedInputStream.newInstance(headerStream);
        codedInput.setSizeLimit(Integer.MAX_VALUE);
        final Alignments.AlignmentHeader header = Alignments.AlignmentHeader.parseFrom(codedInput);


        queryIdentifiers = parseIdentifiers(header.getQueryNameMapping());
        targetIdentifiers = parseIdentifiers(header.getTargetNameMapping());
        if (header.getQueryLengthCount() > 0) {
            queryLengths = new IntArrayList(header.getQueryLengthList()).toIntArray();
        }
        if (header.getTargetLengthCount() > 0) {
            targetLengths = new IntArrayList(header.getTargetLengthList()).toIntArray();
        }
        numberOfQueries = header.getNumberOfQueries();
        numberOfTargets = header.getNumberOfTargets();
        setHeaderLoaded(true);
        numberOfAlignedReads = header.getNumberOfAlignedReads();
    }

    private IndexedIdentifier parseIdentifiers(final Alignments.IdentifierMapping nameMapping) {
        final IndexedIdentifier result = new IndexedIdentifier();

        for (final Alignments.IdentifierInfo info : nameMapping.getMappingsList()) {
            result.put(new MutableString(info.getName()), info.getIndex());
        }
        return result;
    }

    /**
     * {@inheritDoc}
     */
    public void close() {
        alignmentEntryReader.close();
    }

    /**
     * Return the basename corresponding to the input alignment filename.  Note
     * that if the filename does have the extension known to be a compact alignemt
     * the returned value is the original filename
     *
     * @param filename The name of the file to get the basename for
     * @return basename for the alignment file
     */
    public static String getBasename(final String filename) {
        for (final String ext : FileExtensionHelper.COMPACT_ALIGNMENT_FILE_EXTS) {
            if (StringUtils.endsWith(filename, ext)) {
                return StringUtils.removeEnd(filename, ext);
            }
        }

        // perhaps the input was a basename already.
        return filename;
    }

    /**
     * Return the basenames corresponding to the input filenames. Less basename than filenames
     * may be returned (if several filenames reduce to the same baseline after removing
     * the extension).
     *
     * @param filenames The names of the files to get the basnames for
     * @return An array of basenames
     */
    public static String[] getBasenames(final String... filenames) {
        final ObjectSet<String> result = new ObjectArraySet<String>();
        if (filenames != null) {
            for (final String filename : filenames) {
                result.add(getBasename(filename));
            }
        }
        return result.toArray(new String[result.size()]);
    }

    public Iterator<Alignments.AlignmentEntry> iterator() {
        return this;
    }

    public Properties getStatistics() {
        return stats;
    }

    public int getNumberOfAlignedReads() {
        return numberOfAlignedReads;
    }
}
