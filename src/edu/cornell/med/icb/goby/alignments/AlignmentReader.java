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
import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.util.*;
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
    private final FastBufferedMessageChunksReader alignmentEntryReader;
    private Alignments.AlignmentCollection collection;
    private Properties stats;
    private String basename;
    private boolean indexLoaded;
    private int[] targetPositionOffsets;

    /**
     * Required file extension for alignment data in "compact reads" format. Basename + extension
     * must exist and be readable for each extension in this set for an alignment to be
     * readable.
     */
    public static final String[] COMPACT_ALIGNMENT_FILE_REQUIRED_EXTS = {
            ".entries", ".header"
    };

    /**
     * Returns whether this alignment is sorted. Entries in a sorted alignment appear in order of
     * increasing target index and position.
     *
     * @return True if this position is sorted by position. False otherwise.
     */
    public boolean isSorted() {
        return sorted;
    }

    /**
     * Whether this alignment was sorted by position.
     */

    private boolean sorted;

    /**
     * Returns true if this alignment is indexed. When this method returns true,
     * a file called 'basename'.index is expected.
     *
     * @return
     */
    public boolean isIndexed() {
        return indexed;
    }

    private boolean indexed;

    /**
     * Returns true if filename belongs to an alignment basename that can be read.
     *
     * @param filename Filename of an alignment component.
     * @return True if the alignment can be read, false otherwise.
     */
    public static boolean canRead(String filename) {

        String filenameNoExtension = FilenameUtils.removeExtension(filename);
        int count = 0;
        for (String extension : COMPACT_ALIGNMENT_FILE_REQUIRED_EXTS) {
            File fileComponent = new File(filenameNoExtension + extension);

            if (fileComponent.canRead()) {
                // we can read this file.
                count++;
            }
        }
        return count == COMPACT_ALIGNMENT_FILE_REQUIRED_EXTS.length;
    }

    public AlignmentReader(final long startOffset, final long endOffset, final String basename) throws IOException {
        super();
        this.basename = basename;
        final FileInputStream stream = new FileInputStream(basename + ".entries");
        alignmentEntryReader = new FastBufferedMessageChunksReader(startOffset, endOffset, new FastBufferedInputStream(stream));
        LOG.trace("start offset :" + startOffset + " end offset " + endOffset);
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

    public AlignmentReader(final String basename) throws IOException {
        this(0, Long.MAX_VALUE, getBasename(basename));
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

    public AlignmentReader(final InputStream entriesStream) throws IOException {
        super();
        alignmentEntryReader = new FastBufferedMessageChunksReader(0, Long.MAX_VALUE, new FastBufferedInputStream(entriesStream));
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
     * Skip all entries that have position before (targetIndex,position). This method will use the alignment index
     * when it is available to skip directly to the closest chunk start before the entry identified by targetIndex
     * and position.
     *
     * @param targetIndex The index of the target sequence to skip to.
     * @param position    The position on the target sequence.
     * @return The next entry, at position or past position (if not entry at position is found).
     * @throws IOException If an error occurs reading the alignment header. The header is accessed to check that the alignment is sorted.
     */
    public final Alignments.AlignmentEntry skipTo(final int targetIndex, final int position) throws IOException {

        reposition(targetIndex, position);
        Alignments.AlignmentEntry entry = null;
        boolean hasNext = false;
        while ((hasNext = hasNext()) &&
                (entry = next()) != null &&
                (entry.getTargetIndex() < targetIndex ||
                        (entry.getTargetIndex() == targetIndex && entry.getPosition() < position))) {
        }

        if (!hasNext) return null;
        else return entry;

    }

    /**
     * Reposition the reader to a new target sequence and start position.
     *
     * @param targetIndex Index of the target sequence to reposition to.
     * @param position    Position in the target sequence to reposition to.
     * @throws IOException If an error occurs repositioning.
     */
    public final void reposition(final int targetIndex, final int position) throws IOException {
        readHeader();
        if (!sorted) throw new UnsupportedOperationException("skipTo cannot be used with unsorted alignments.");

        readIndex();
        repositionInternal(targetIndex, position);
    }

    private void repositionInternal(final int targetIndex, final int position) throws IOException {
        if (!indexLoaded) return;
        final int absolutePosition = recodePosition(targetIndex, position);
        int offsetIndex = Arrays.binarySearch(indexAbsolutePositions.elements(), absolutePosition);
        offsetIndex = offsetIndex < 0 ? -1 - offsetIndex : offsetIndex;
        offsetIndex = offsetIndex >= indexOffsets.size() ? indexOffsets.size() - 1 : offsetIndex;
        if (offsetIndex < 0) {
            // empty alignment.
            return;
        }

        final long newPosition = indexOffsets.getLong(offsetIndex);
        final long currentPosition = alignmentEntryReader.position();
        if (newPosition > currentPosition) {

            alignmentEntryReader.seek(newPosition);
        }
    }

    protected int recodePosition(final int firstTargetIndexInChunk, final int firstPositionInChunk) {
        return targetPositionOffsets[firstTargetIndexInChunk] + firstPositionInChunk;
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
        if (!isHeaderLoaded()) {
            // accept very large header messages, since these may contain query identifiers:
            final CodedInputStream codedInput = CodedInputStream.newInstance(headerStream);
            codedInput.setSizeLimit(Integer.MAX_VALUE);
            final Alignments.AlignmentHeader header = Alignments.AlignmentHeader.parseFrom(codedInput);

            smallestQueryIndex = header.getSmallestSplitQueryIndex();
            largestQueryIndex = header.getLargestSplitQueryIndex();
            queryIdentifiers = parseIdentifiers(header.getQueryNameMapping());
            targetIdentifiers = parseIdentifiers(header.getTargetNameMapping());
            if (!header.getQueryLengthsStoredInEntries()) {
                if (header.hasConstantQueryLength()) {
                    this.constantQueryLengths = true;
                    this.constantLength = header.getConstantQueryLength();
                    queryLengths = null;
                } else if (header.getQueryLengthCount() > 0) {
                    queryLengths = new IntArrayList(header.getQueryLengthList()).toIntArray();
                    compactQueryLengths();
                }
            } else {
                queryLengths = null;
                constantQueryLengths = false;
            }
            if (header.getTargetLengthCount() > 0) {
                targetLengths = new IntArrayList(header.getTargetLengthList()).toIntArray();
            }
            numberOfQueries = header.getNumberOfQueries();
            numberOfTargets = header.getNumberOfTargets();
            setHeaderLoaded(true);
            numberOfAlignedReads = header.getNumberOfAlignedReads();
            sorted = header.getSorted();
            indexed = header.getIndexed();

        }
    }

    private LongArrayList indexOffsets = new LongArrayList();
    private LongArrayList indexAbsolutePositions = new LongArrayList();

    /**
     * Read the index. The header is also loaded.
     *
     * @throws IOException If an error occurs loading the index or header.
     */
    public final void readIndex() throws IOException {
        if (indexed && !indexLoaded) {
            // header is needed to access target lengths:
            readHeader();
            GZIPInputStream indexStream = new GZIPInputStream(new FileInputStream(basename + ".index"));

            final CodedInputStream codedInput = CodedInputStream.newInstance(indexStream);
            codedInput.setSizeLimit(Integer.MAX_VALUE);
            final Alignments.AlignmentIndex index = Alignments.AlignmentIndex.parseFrom(codedInput);
            indexOffsets.clear();
            indexAbsolutePositions.clear();

            for (long offset : index.getOffsetsList()) {
                indexOffsets.add(offset);
            }
            for (long absolutePosition : index.getAbsolutePositionsList()) {
                indexAbsolutePositions.add(absolutePosition);
            }
            // trimming is essential for the binary search to work reliably with the result of elements():
            indexAbsolutePositions.trim();
            indexOffsets.trim();
// calculate the coding offset for each target index. This information will be used by recode
            targetPositionOffsets = new int[targetLengths.length];
            for (int targetIndex = 0; targetIndex < targetLengths.length; targetIndex++) {
                targetPositionOffsets[targetIndex] += targetLengths[targetIndex];
                targetPositionOffsets[targetIndex] += targetIndex < 1 ? 0 : targetPositionOffsets[targetIndex - 1];
            }

            indexLoaded = true;
        }
    }

    private void compactQueryLengths() {
        final IntSet distinctQueryLength = new IntOpenHashSet();

        for (final int length : queryLengths) {
            if (length != 0) {
                distinctQueryLength.add(length);
            }
        }
        if (distinctQueryLength.size() == 1) {
            this.constantQueryLengths = true;
            this.constantLength = distinctQueryLength.iterator().nextInt();
            //release the space:
            queryLengths = null;
        }
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

    /**
     * Returns whether this read has query length information.
     *
     * @return True or false.
     */
    public boolean hasQueryLengths() {
        return ArrayUtils.isNotEmpty(queryLengths);
    }
}
