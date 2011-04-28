/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.ArrayUtils;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * Reads alignments written with AlignmentWriter.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 6:36:04 PM
 */
public class AlignmentReaderImpl extends AbstractAlignmentReader implements AlignmentReader {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(AlignmentReaderImpl.class);

    private InputStream headerStream;
    private int numberOfAlignedReads;
    private final FastBufferedMessageChunksReader alignmentEntryReader;
    private Alignments.AlignmentCollection collection;
    private Properties stats;
    private String basename;
    private boolean indexLoaded;
    private long[] targetPositionOffsets;

    private int endReferenceIndex;
    private int endPosition;
    private int startPosition;
    private int startReferenceIndex;
    /**
     * Required file extension for alignment data in "compact reads" format. Basename + extension
     * must exist and be readable for each extension in this set for an alignment to be
     * readable.
     */
    public static final String[] COMPACT_ALIGNMENT_FILE_REQUIRED_EXTS = {
            ".entries", ".header"
    };
    private Alignments.AlignmentEntry nextEntry;
    private Alignments.AlignmentEntry nextEntryNoFilter;
    private boolean queryLengthStoredInEntries;
    private String alignerName;
    private String alignerVersion;


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
    public static boolean canRead(final String filename) {

        final String filenameNoExtension = FilenameUtils.removeExtension(filename);
        int count = 0;
        for (final String extension : AlignmentReaderImpl.COMPACT_ALIGNMENT_FILE_REQUIRED_EXTS) {
            final File fileComponent = new File(filenameNoExtension + extension);

            if (fileComponent.canRead()) {
                // we can read this file.
                count++;
            }
        }
        return count == AlignmentReaderImpl.COMPACT_ALIGNMENT_FILE_REQUIRED_EXTS.length;
    }

    /**
     * A constructor that allows reading a slice of an alignment file contained exactly between a start
     * and an end location. Start and end locations are genomic/reference positions. Entries will be returned
     * that occur after the start position and before the end position.
     *
     * @param basename            Basename for the alignemnt.
     * @param startReferenceIndex Index of the reference for the start position.
     * @param startPosition       Position on the reference for the start position.
     * @param endReferenceIndex   Index of the reference for the end position.
     * @param endPosition         Position on the reference for the end position.
     * @throws IOException Thrown if an error occurs opening or reading the alignment file.
     */
    public AlignmentReaderImpl(final String basename,
                               final int startReferenceIndex,
                               final int startPosition,
                               final int endReferenceIndex,
                               final int endPosition)
            throws IOException {

        super();
        this.basename = basename;

        try {
            headerStream = new GZIPInputStream(new FileInputStream(basename + ".header"));
        } catch (IOException e) {
            // try not compressed for compatibility with 1.4-:
            LOG.trace("falling back to legacy 1.4- uncompressed header.");

            headerStream = new FileInputStream(basename + ".header");
        }

        readHeader();
        if (!indexed)
            throw new UnsupportedOperationException("The alignment must be sorted and indexed to read slices of data by reference position.");
        readIndex();
        final FileInputStream stream = new FileInputStream(basename + ".entries");
        final long startOffset = getByteOffset(startReferenceIndex, startPosition, 0);
        long endOffset = getByteOffset(endReferenceIndex, endPosition + 1, 1);


        this.endPosition = endPosition;
        this.endReferenceIndex = endReferenceIndex;
        this.startPosition = startPosition;
        this.startReferenceIndex = startReferenceIndex;
        alignmentEntryReader = new FastBufferedMessageChunksReader(startOffset > 0 ? startOffset : 0,
                endOffset > 0 ? endOffset : Long.MAX_VALUE,
                new FastBufferedInputStream(stream));
        LOG.trace("start offset :" + startOffset + " end offset " + endOffset);

        stats = new Properties();
        final File statsFile = new File(basename + ".stats");
        if (statsFile.exists()) {
            Reader statsFileReader = null;
            try {
                statsFileReader = new FileReader(statsFile);
                stats.load(statsFileReader);
            } catch (IOException e) {
                LOG.warn("cannot load properties for basename: " + basename, e);
                throw e;
            } finally {
                IOUtils.closeQuietly(statsFileReader);
            }
        }
    }


    /**
     * Open a Goby alignment file for reading between the byte positions startOffset and endOffset.
     *
     * @param startOffset Position in the file where reading will start (in bytes).
     * @param endOffset   Position in the file where reading will end (in bytes).
     * @param basename    Basename of the alignment to read.
     * @throws IOException If an error occurs opening or reading the file.
     */
    public AlignmentReaderImpl(final long startOffset, final long endOffset, final String basename) throws IOException {
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
        startReferenceIndex = 0;
        startPosition = 0;
        endReferenceIndex = Integer.MAX_VALUE;
        endPosition = Integer.MAX_VALUE;
    }

    public AlignmentReaderImpl(final String basename) throws IOException {
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

    public AlignmentReaderImpl(final InputStream entriesStream) throws IOException {
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
    public AlignmentReaderImpl(final long start, final long end, final FastBufferedInputStream stream)
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
        //   System.out.println("hasNext");
        if (nextEntry != null) return true;

        int entryTargetIndex;
        int position;
        do {

            if (!hasNextEntry()) {
                nextEntry = null;
                return false;
            }

            nextEntry = nextEntry();
            entryTargetIndex = nextEntry.getTargetIndex();
            // Early stop if we are past the end location:
            position = nextEntry.getPosition();
            if (entryTargetIndex > endReferenceIndex ||
                    (entryTargetIndex == endReferenceIndex && position > endPosition)) {
                nextEntry = null;
                return false;
            }
        } while (entryTargetIndex < startReferenceIndex ||
                (entryTargetIndex == startReferenceIndex && position < startPosition));


        return true;
    }

    /**
     * Returns the next alignment entry from the input stream.
     *
     * @return the alignment read entry from the input stream.
     */
    public Alignments.AlignmentEntry next() {
        // System.out.println("next");
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        try {
            if (LOG.isTraceEnabled()) {
                LOG.trace(String.format("Returning next entry at position %s/%d%n", getTargetIdentifiers().get(new MutableString(nextEntry.getTargetIndex())), nextEntry.getPosition()));
            }
            return nextEntry;

        } finally {

            nextEntry = null;
        }

    }

    /**
     * Returns true if the input has more entries.
     *
     * @return true if the input has more entries, false otherwise.
     */
    private boolean hasNextEntry() {
        //    System.out.println("hasNextEntry");
        if (nextEntryNoFilter != null) return true;

        if (collection != null && alignmentEntryReader.getEntryIndex() < collection.getAlignmentEntriesCount()) {
            nextEntryNoFilter = collection.getAlignmentEntries(alignmentEntryReader.getEntryIndex());
            alignmentEntryReader.incrementEntryIndex();
            return true;
        } else {
            collection = null;
            final boolean hasNext = alignmentEntryReader.hasNext(collection, numberOfEntries());

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
            if (hasNext) {
                nextEntryNoFilter = collection.getAlignmentEntries(alignmentEntryReader.getEntryIndex());
                alignmentEntryReader.incrementEntryIndex();
                return true;
            } else return false;
        }


    }

    private Alignments.AlignmentEntry nextEntry() {
        //      System.out.println("nextEntry");
        if (!hasNextEntry()) {
            throw new NoSuchElementException();
        }
        try {
            return nextEntryNoFilter;
        } finally {
            nextEntryNoFilter = null;
        }
        //    LOG.warn(String.format("returning entry %d/%d%n", entry.getTargetIndex(), entry.getPosition()));
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
        /*
        If the skipTo requested position is out of the position window provided to the constructor, return
        no result. Otherwise, adjust the position to be within the window.
         */
        if (targetIndex < this.startReferenceIndex || targetIndex > this.endReferenceIndex) {
            return null;
        }
        int positionChanged = position;
        if (targetIndex == this.startReferenceIndex) {
            positionChanged = Math.max(this.startPosition, position);
        }
        if (LOG.isTraceEnabled()) {
            LOG.trace(String.format("skipTo %d/%d%n", targetIndex, positionChanged));
        }
        reposition(targetIndex, positionChanged);
        Alignments.AlignmentEntry entry = null;
        boolean hasNext = false;
        while ((hasNext = hasNext()) &&
                ((entry = next()) != null) &&
                (entry.getTargetIndex() < targetIndex ||
                        (entry.getTargetIndex() == targetIndex && entry.getPosition() < positionChanged))) {
        }

        if (!hasNext) {
            return null;
        } else {
            return entry;
        }
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
        if (!sorted) {
            throw new UnsupportedOperationException("skipTo cannot be used with unsorted alignments.");
        }

        readIndex();
        repositionInternal(targetIndex, position);
    }

    private void repositionInternal(final int targetIndex, final int position) throws IOException {
        if (!indexLoaded) {
            return;
        }
        final long absolutePosition = recodePosition(targetIndex, position);
        int offsetIndex = Arrays.binarySearch(indexAbsolutePositions.elements(), absolutePosition);
        offsetIndex = offsetIndex < 0 ? -1 - offsetIndex : offsetIndex;
        // NB offsetIndex contains absolutePosition in the first entry, but the chunk before it also likely
        // contains entries with this absolute position. We therefore substract one to position on the chunk
        // before. 
        offsetIndex = offsetIndex >= indexOffsets.size() ? indexOffsets.size() - 1 : offsetIndex - 1;

        if (offsetIndex < 0) {
            // empty alignment.
            return;
        }
        // if (indexAbsolutePositions.getLong(offsetIndex)<absolutePosition) {
        // the offset was not found, indicating that the chunk starts before (targetIndex,position) location.
        // This means that we must scan the chunk immediately before to check for entries that can be at this
        // position as well.
        //    --offsetIndex;
        //   }
        final long newPosition = indexOffsets.getLong(offsetIndex);
        final long currentPosition = alignmentEntryReader.position();
        if (newPosition > currentPosition) {

            alignmentEntryReader.seek(newPosition);
        }
    }

    /**
     * Calculate the offset (in bytes) in the compact entries file for a specific targetIndex and position.
     * Entries that can be read after this position are garanteed to have targetIndex larger or equal to targetIndex
     * and positions larger or equal to position.
     * The parameter chunkOffset can be used to iterate through successive protocol buffer compressed chunks.
     * A typical usage is to call getByteOffset with startReference and startPosition. A second call to getByteOffset
     * with endReference and endPosition with chunkOffset=0 will return an end position. If the end and start positions
     * are the same, both start and end locations are in the same chunk. In this case, it is necessary to extend the end
     * byte position to include the entire chunk. This can be achieved by calling getByteOffset with  endReference and endPosition
     * and a chunkOffset of 1. Because it is possible that the next chunk also contains only entries with the same position,
     * one must
     * use 0, 1,2, etc until the returned offset differs allows to make sure the indexed position is at the start of a new chunk.
     *
     * @param targetIndex Index of a reference sequence.
     * @param position    Position along the reference sequence.
     * @param chunkOffset Offset used to iterate through successive PB chunks.
     * @return the largest position in byte in the entries file that occur before the location (targetIndex, position) or
     *         Long.MIN_VALUE if the offset cannot be determined (e.g., alignment is empty).
     */
    protected long getByteOffset(final int targetIndex, final int position, final int chunkOffset) {

        if (targetIndex > targetPositionOffsets.length) return Long.MAX_VALUE;

        final long absolutePosition = recodePosition(targetIndex, position);
        int offsetIndex = Arrays.binarySearch(indexAbsolutePositions.elements(), absolutePosition);
        offsetIndex = offsetIndex < 0 ? -1 - offsetIndex : offsetIndex;
        offsetIndex = offsetIndex >= indexOffsets.size() ? indexOffsets.size() - 1 : offsetIndex - 1;
        if (offsetIndex < 0) {
            // empty alignment.
            return Long.MIN_VALUE;
        }

        if (offsetIndex + chunkOffset < indexOffsets.size()) {
            final long byteOffset = indexOffsets.getLong(offsetIndex + chunkOffset);
            return byteOffset;
        } else {
            // return an end-offset past the beginning of the last chunk:
            return indexOffsets.getLong(offsetIndex) + 10;
        }

    }

    protected long recodePosition(final int firstTargetIndexInChunk, final int firstPositionInChunk) {
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

            alignerName = header.getAlignerName();
            alignerVersion = header.getAlignerVersion();
            smallestQueryIndex = header.getSmallestSplitQueryIndex();
            largestQueryIndex = header.getLargestSplitQueryIndex();
            queryIdentifiers = parseIdentifiers(header.getQueryNameMapping());
            targetIdentifiers = parseIdentifiers(header.getTargetNameMapping());
            if (header.hasConstantQueryLength()) {
                this.constantQueryLengths = true;
                this.constantLength = header.getConstantQueryLength();
            }
            queryLengthStoredInEntries = header.getQueryLengthsStoredInEntries();

            assert queryLengthStoredInEntries : "This version of Goby requires that query lengths are stored in entries." +
                    " You can upgrade old alignment files by transfering data with the concat mode of a previous version.";

            if (header.getTargetLengthCount() > 0) {
                targetLengths = new IntArrayList(header.getTargetLengthList()).toIntArray();
            }
            numberOfQueries = header.getNumberOfQueries();
            numberOfTargets = header.getNumberOfTargets();
            setHeaderLoaded(true);
            numberOfAlignedReads = header.getNumberOfAlignedReads();
            // we determine sortedness and index state from the header and by checking that the index file exists.
            // This allows to recover alignments when the index file was deleted. We can then read and sort them
            // again.
            sorted = header.getSorted() && indexExists(basename);
            indexed = header.getIndexed()&& indexExists(basename);

        }
    }

    private boolean indexExists(String basename) {
        return new File(basename+".index").exists();
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
            final GZIPInputStream indexStream = new GZIPInputStream(new FileInputStream(basename + ".index"));

            final CodedInputStream codedInput = CodedInputStream.newInstance(indexStream);
            codedInput.setSizeLimit(Integer.MAX_VALUE);
            final Alignments.AlignmentIndex index = Alignments.AlignmentIndex.parseFrom(codedInput);
            indexOffsets.clear();
            indexAbsolutePositions.clear();

            for (final long offset : index.getOffsetsList()) {
                indexOffsets.add(offset);
            }
            for (final long absolutePosition : index.getAbsolutePositionsList()) {
                indexAbsolutePositions.add(absolutePosition);
            }
            // trimming is essential for the binary search to work reliably with the result of elements():
            indexAbsolutePositions.trim();
            indexOffsets.trim();
// calculate the coding offset for each target index. This information will be used by recode
            targetPositionOffsets = new long[targetLengths.length];
            for (int targetIndex = 0; targetIndex < targetLengths.length; targetIndex++) {
                targetPositionOffsets[targetIndex] += targetLengths[targetIndex];
                targetPositionOffsets[targetIndex] += targetIndex < 1 ? 0 : targetPositionOffsets[targetIndex - 1];
            }

            indexLoaded = true;
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
        if (alignmentEntryReader != null) {
            alignmentEntryReader.close();
        }
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
     * the extension). This method returns the unique set of basenames in the same order they are
     * provided as argument.
     *
     * @param filenames The names of the files to get the basnames for
     * @return An array of basenames
     */
    public static String[] getBasenames(final String... filenames) {
        final ObjectSet<String> result = new ObjectArraySet<String>();
        final ObjectList<String> unique = new ObjectArrayList<String>();
        if (filenames != null) {
            for (final String filename : filenames) {


                final String newBasename = getBasename(filename);
                if (!result.contains(newBasename)) {
                    unique.add(newBasename);
                    result.add(getBasename(filename));
                }
            }
        }
        return unique.toArray(new String[unique.size()]);
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
     * Returns a sample of locations covered by this alignment.
     *
     * @param modulo Modulo to avoid sampling every position in the genome.
     * @return A set of positions that do occur in the genome, rounded to the specified modulo value (absoluteLocation-(absoluteLocation % modulo)).
     *         * @throws IOException
     */
    public ObjectList<ReferenceLocation> getLocations(int modulo) throws IOException {
        if (!isIndexed()) throw new RuntimeException("Alignment must be sorted and indexed to obtain locations.");
        readIndex();
        ObjectList<ReferenceLocation> result = new ObjectArrayList<ReferenceLocation>();
        for (long absoluteLocation : indexAbsolutePositions) {
            final ReferenceLocation location = decodeAbsoluteLocation(absoluteLocation - (absoluteLocation % modulo));
            result.add(location);

        }
        return result;
    }

    private ReferenceLocation decodeAbsoluteLocation(long absoluteLocation) {
        int referenceIndex;
        for (referenceIndex = this.targetPositionOffsets.length - 1; referenceIndex >= 0; referenceIndex--) {
            final long offset = this.targetPositionOffsets[referenceIndex];
            if (absoluteLocation > offset) {
                absoluteLocation -= offset;
                break;
            }
        }

        return new ReferenceLocation(referenceIndex + 1, (int) absoluteLocation);
    }

    public boolean isQueryLengthStoredInEntries() {
        return queryLengthStoredInEntries;
    }

    /**
     * Return the name of the aligner that produced this alignment.
     *
     * @return the name of the aligner that produced this alignment.
     */
    public String getAlignerName() {
        return alignerName;
    }

    /**
     * Return the version of the aligner that produced this alignment.
     *
     * @return the version of the aligner that produced this alignment.
     */
    public String getAlignerVersion() {
        return alignerVersion;
    }

    public int getConstantQueryLength() {
        return constantLength;
    }
}
