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
import edu.cornell.med.icb.goby.compression.ChunkCodec;
import edu.cornell.med.icb.goby.compression.FastBufferedMessageChunksReader;
import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

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
    /**
     * Other possible extensions that can follow a Goby alignment basename.
     */
    public static final String[] COMPACT_ALIGNMENT_FILE_POSSIBLE_EXTS = {
            ".index", ".perm", ".tmh"
    };
    private Alignments.AlignmentEntry nextEntry;
    private Alignments.AlignmentEntry nextEntryNoFilter;
    private boolean queryLengthStoredInEntries;
    private String alignerName;
    private String alignerVersion;
    /**
     * The version of Goby that created the alignment file we are reading.
     */
    private String gobyVersion;
    private boolean queryIndicesWerePermuted;
    /**
     * This field is true when all the entries of the reader have the read_quality_scores field populated.
     */
    private boolean allReadQualityScores;
    /**
     * Start offset of the slice this read was created with. Number of bytes into the entries file were to start reading from.
     */
    private long startOffset;
    /**
     * End offset of the slice this read was created with. Number of bytes into the entries file were to stop reading.
     */
    private long endOffset;
    private boolean hasQueryIndexOccurrences;
    private List<Alignments.ReadOriginInfo> readOriginInfoList;
    private boolean hasAmbiguity;


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
        String fileExtension = FilenameUtils.getExtension(filename);

        if (!(ArrayUtils.contains(AlignmentReaderImpl.COMPACT_ALIGNMENT_FILE_POSSIBLE_EXTS, "." + fileExtension) ||
                ArrayUtils.contains(AlignmentReaderImpl.COMPACT_ALIGNMENT_FILE_REQUIRED_EXTS, "." + fileExtension))) {
            // the file does not contain any of the possible Goby extensions. It is not a supported file.
            return false;
        }
        // the file contains a Goby alignment extension, we further check that each needed extension exists:
        int count = 0;
        for (final String extension : AlignmentReaderImpl.COMPACT_ALIGNMENT_FILE_REQUIRED_EXTS) {

            if (RepositionableInputStream.resourceExist(filenameNoExtension + extension)) {
                // we can read this file.
                count++;
            }
        }
        return count == AlignmentReaderImpl.COMPACT_ALIGNMENT_FILE_REQUIRED_EXTS.length;
    }

    /**
     * A constructor that allows reading a slice of an alignment file contained exactly between a start
     * and an end location. Start and end locations are genomic/reference positions. Entries will be returned
     * that occur after the start position and up to the end position (start and end positions are inclusive).
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
        this(basename,
                startReferenceIndex,
                startPosition,
                endReferenceIndex,
                endPosition, false);
    }

    /**
     * A constructor that allows reading a slice of an alignment file contained exactly between a start
     * and an end location. Start and end locations are genomic/reference positions. Entries will be returned
     * that occur after the start position and up to the end position (start and end positions are inclusive).
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
                               final int endPosition, final boolean upgrade)
            throws IOException {

        super(true, getBasename(basename));
        this.basename = getBasename(basename);

        try {
            headerStream = new GZIPInputStream(new RepositionableInputStream((this.basename + ".header")));
        } catch (IOException e) {
            // try not compressed for compatibility with 1.4-:
            LOG.trace("falling back to legacy 1.4- uncompressed header.");

            headerStream = new RepositionableInputStream(this.basename + ".header");
        }

        readHeader();
        if (!indexed)
            throw new UnsupportedOperationException("The alignment must be sorted and indexed to read slices of data by reference position.");
        readIndex();
        final InputStream stream = new RepositionableInputStream(this.basename + ".entries");
        final long startOffset = getByteOffset(startReferenceIndex, startPosition, 0);
        long endOffset = getByteOffset(endReferenceIndex, endPosition + 1, 1);


        this.endPosition = endPosition;
        this.endReferenceIndex = endReferenceIndex;
        this.startPosition = startPosition;
        this.startReferenceIndex = startReferenceIndex;
        alignmentEntryReader = new FastBufferedMessageChunksReader(startOffset > 0 ? startOffset : 0,
                endOffset > 0 ? endOffset : Long.MAX_VALUE,
                new FastBufferedInputStream(stream));
        alignmentEntryReader.setHandler(new AlignmentCollectionHandler());
        LOG.trace("start offset :" + startOffset + " end offset " + endOffset);

        stats = new Properties();
        String statsFilename = basename + ".stats";

        if (RepositionableInputStream.resourceExist(statsFilename)) {
            final File statsFile = new File(statsFilename);
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
     * Please note that this method does not attempt to upgrade the alignment since the offsets would likely be wrong
     * against an upgraded file.
     *
     * @param startOffset Position in the file where reading will start (in bytes).
     * @param endOffset   Position in the file where reading will end (in bytes).
     * @param basename    Basename of the alignment to read.
     * @throws IOException If an error occurs opening or reading the file.
     */
    public AlignmentReaderImpl(final long startOffset, final long endOffset, final String basename) throws IOException {
        // do not try to upgrade when we provide explicit offsets since the offsets would be wrong in the upgraded file anyway.
        this(startOffset, endOffset, basename, false);
    }

    /**
     * Open a Goby alignment file for reading between the byte positions startOffset and endOffset.
     *
     * @param startOffset Position in the file where reading will start (in bytes).
     * @param endOffset   Position in the file where reading will end (in bytes).
     * @param basename    Basename of the alignment to read.
     * @param upgrade     Whether to try to upgrade the alignment on the fly.
     * @throws IOException If an error occurs opening or reading the file.
     */
    public AlignmentReaderImpl(final long startOffset, final long endOffset, final String basename, boolean upgrade) throws IOException {
        super(upgrade, getBasename(basename));
        this.basename = getBasename(basename);
        this.startOffset = startOffset;
        this.endOffset = endOffset;
        final String entriesFile = this.basename + ".entries";
        boolean entriesFileExist = RepositionableInputStream.resourceExist(entriesFile);
        if (entriesFileExist) {
            final InputStream stream = new RepositionableInputStream(entriesFile);

            alignmentEntryReader = new FastBufferedMessageChunksReader(startOffset, endOffset, new FastBufferedInputStream(stream));
            alignmentEntryReader.setHandler(new AlignmentCollectionHandler());
        } else {
            alignmentEntryReader = null;
        }
        LOG.trace("start offset :" + startOffset + " end offset " + endOffset);
        try {
            headerStream = new GZIPInputStream(new RepositionableInputStream(this.basename + ".header"));
        } catch (IOException e) {
            e.printStackTrace();
            // try not compressed for compatibility with 1.4-:
            LOG.trace("falling back to legacy 1.4- uncompressed header.");

            headerStream = new RepositionableInputStream(this.basename + ".header");
        }
        stats = new Properties();
        String statsFilename = this.basename + ".stats";

        if (RepositionableInputStream.resourceExist(statsFilename)) {
            final File statsFile = new File(statsFilename);
            Reader statsFileReader = null;
            try {
                statsFileReader = new InputStreamReader(new RepositionableInputStream(statsFilename));
                stats.load(statsFileReader);
            } catch (IOException e) {
                LOG.warn("cannot load properties for basename: " + this.basename, e);
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
        this(0, Long.MAX_VALUE, getBasename(basename), true);

    }

    public AlignmentReaderImpl(final String basename, boolean upgrade) throws IOException {
        this(0, Long.MAX_VALUE, getBasename(basename), upgrade);

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
        super(true, null);
        alignmentEntryReader = new FastBufferedMessageChunksReader(0, Long.MAX_VALUE, new FastBufferedInputStream(entriesStream));
        alignmentEntryReader.setHandler(new AlignmentCollectionHandler());
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
        super(true, null);
        alignmentEntryReader = new FastBufferedMessageChunksReader(start, end, stream);
        alignmentEntryReader.setHandler(new AlignmentCollectionHandler());
    }

    private int numberOfEntries() {
        return collection != null ? collection.getAlignmentEntriesCount() : 0;
    }

    private IndexedIdentifier identifiers;
    private DoubleIndexedIdentifier back;

    /**
     * Returns true if the input has more entries.
     *
     * @return true if the input has more entries, false otherwise.
     */
    public boolean hasNext() {

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
                nextEntryNoFilter = null;
                collection = null;
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

            if (!nextEntry.hasMultiplicity()) {
                // set the default multiplicity when the field was not defined.
                Alignments.AlignmentEntry.Builder builder = Alignments.AlignmentEntry.newBuilder(nextEntry);
                builder.setMultiplicity(1);
                nextEntry = builder.build();
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

            final ChunkCodec codec = alignmentEntryReader.getChunkCodec();
            try {
                final byte[] compressedBytes = alignmentEntryReader.getCompressedBytes();
                if (compressedBytes != null) {
                    collection = (Alignments.AlignmentCollection) codec.decode(compressedBytes);
                    if (collection == null || collection.getAlignmentEntriesCount() == 0) {
                        return false;
                    }
                    if (LOG.isTraceEnabled()) {
                        if (back == null) {
                            readHeader();
                            identifiers = getTargetIdentifiers();
                            back = new DoubleIndexedIdentifier(identifiers);
                        }

                        final Alignments.AlignmentEntry firstEntry = collection.getAlignmentEntries(0);
                        if (targetPositionOffsets != null) {
                            LOG.trace(String.format("New collection with first entry at position id=%s/pos=%d absolutePosition=%d %n", back.getId(firstEntry.getTargetIndex()),
                                    firstEntry.getPosition(),
                                    recodePosition(firstEntry.getTargetIndex(), firstEntry.getPosition())
                            ));
                        }

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
            readHeader();
            if (back == null) {

                identifiers = getTargetIdentifiers();
                back = new DoubleIndexedIdentifier(identifiers);
            }
            if (targetPositionOffsets != null) {
                LOG.trace(String.format("skipTo id=%s/pos=%d absolutePosition=%d %n", back.getId(targetIndex), positionChanged,
                        recodePosition(targetIndex, position)
                ));
            }

        }

        repositionInternal(targetIndex, positionChanged, false);
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
            throw new UnsupportedOperationException("reposition cannot be used with unsorted alignments.");
        }
        readIndex();
        this.alignmentEntryReader.flush();
        repositionInternal(targetIndex, position, true);
    }

    /**
     * Reposition to a genomic position. The goBack flag, when true, allows to reposition to positions that
     * we have already passed. When false, reposition will only advance to future positions.
     *
     * @param targetIndex index of the target sequence to reposition to.
     * @param position    position of the location to reposition to.
     * @param goBack      If true, the method will reposition to past positions, otherwise, only reposition to future locations.
     * @throws IOException If an error occurs repositioning.
     */
    private void repositionInternal(final int targetIndex, final int position, boolean goBack) throws IOException {
//        assert indexLoaded : "index must be loaded to repositionInternal.";
        if (!indexLoaded) {
            return;
        }
        final long absolutePosition = recodePosition(targetIndex, position);
        int offsetIndex = Arrays.binarySearch(indexAbsolutePositions.elements(), absolutePosition);
        offsetIndex = offsetIndex < 0 ? -1 - offsetIndex : offsetIndex;
        // NB offsetIndex contains absolutePosition in the first entry, but the chunk before it also likely
        // contains entries with this absolute position. We therefore substract one to position on the chunk
        // before.
        offsetIndex = offsetIndex >= indexOffsets.size() ? indexOffsets.size() - 1 : Math.max(offsetIndex - 1, 0);


        if (offsetIndex < 0) {
            // empty alignment.
            return;
        }

        // max below ensures we never go back to before the start of the slice the reader was restricted to at
        // construction time:
        final long newBytePosition = Math.max(startOffset, indexOffsets.getLong(offsetIndex));
        final long currentPosition = alignmentEntryReader.position();
        if (newBytePosition >= currentPosition) {

            seek(newBytePosition);
        } else {
            if (goBack) {
                // only reposition to past locations if we are called directly through reposition. Otherwise, we honor
                // the skipTo contract and do not reposition to previously visited locations.
                seek(newBytePosition);
            }
        }
    }

    /**
     * Seek to a position in the entries file.
     *
     * @param byteOffset Where to reposition the stream.
     * @throws IOException If an error occured.
     */
    protected void seek(long byteOffset) throws IOException {
        alignmentEntryReader.seek(byteOffset);
        nextEntry = null;
        nextEntryNoFilter = null;
        collection = null;
    }

    /**
     * Calculate the offset (in bytes) in the compact entries file for a specific targetIndex and position.
     * Entries that can be read after this position are guaranteed to have targetIndex larger or equal to targetIndex
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

        if (targetIndex >= targetPositionOffsets.length) return Long.MAX_VALUE;

        final long absolutePosition = recodePosition(targetIndex, position);
        int offsetIndex = Arrays.binarySearch(indexAbsolutePositions.elements(), absolutePosition);
        offsetIndex = offsetIndex < 0 ? -1 - offsetIndex : offsetIndex;
        offsetIndex = offsetIndex >= indexOffsets.size() ? indexOffsets.size() - 1 : offsetIndex - 1;
        if (offsetIndex + chunkOffset < 0) {
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
            queryIndicesWerePermuted = header.getQueryIndicesWerePermuted();

            if (header.getTargetLengthCount() > 0) {
                targetLengths = new IntArrayList(header.getTargetLengthList()).toIntArray();
            }
            numberOfQueries = header.getNumberOfQueries();
            numberOfTargets = header.getNumberOfTargets();

            numberOfAlignedReads = header.getNumberOfAlignedReads();
            // we determine sortedness and index state from the header and by checking that the index file exists.
            // This allows to recover alignments when the index file was deleted. We can then read and sort them
            // again.
            sorted = header.getSorted() && indexExists(basename);
            indexed = header.getIndexed() && indexExists(basename);
            gobyVersion = header.getVersion();
            allReadQualityScores = header.getAllReadQualityScores();
            hasQueryIndexOccurrences = header.getQueryIndexOccurrences();
            readOriginInfoList = header.getReadOriginList();
            hasAmbiguity = header.getAmbiguityStoredInEntries();

            setHeaderLoaded(true);
        }
    }

    /**
     * Indicates if all the entries of this alignment have the read_quality_score field.
     *
     * @return True or False.
     */
    public boolean getHasAllReadQualityScores() {
        return allReadQualityScores;
    }

    private boolean indexExists(String basename) {
        return RepositionableInputStream.resourceExist(basename + ".index");
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
            final GZIPInputStream indexStream = new GZIPInputStream(new RepositionableInputStream(basename + ".index"));

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
            targetPositionOffsets[0] = 0;
            for (int targetIndex = 1; targetIndex < targetLengths.length; targetIndex++) {
                targetPositionOffsets[targetIndex] =
                        targetLengths[targetIndex - 1] +
                                targetPositionOffsets[targetIndex - 1];

            }
            indexLoaded = true;
        }
        if (!indexLoaded && !indexed) {
            LOG.warn("Trying to read index for an alignment that is not indexed.");
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


    public Iterator<Alignments.AlignmentEntry> iterator() {
        return this;
    }

    public Properties getStatistics() {
        return stats;
    }

    public int getNumberOfAlignedReads() {
        return numberOfAlignedReads;
    }


    @Override
    public ObjectList<ReferenceLocation> getLocationsByBytes(int bytesPerSlice) throws IOException {

        assert isHeaderLoaded() : "header must be loaded to query locations.";
        if (!isIndexed()) throw new RuntimeException("Alignment must be sorted and indexed to obtain locations.");

        readIndex();
        int i = 0;
        ObjectList<ReferenceLocation> result = new ObjectArrayList<ReferenceLocation>();
        long lastFileOffsetPushed = -1;

        for (long absoluteLocation : indexAbsolutePositions) {
            long offsetInEntriesFile = indexOffsets.get(i);
            long compressedByteAmountSincePreviousLocation = offsetInEntriesFile - lastFileOffsetPushed;
            if (lastFileOffsetPushed==-1 || compressedByteAmountSincePreviousLocation > bytesPerSlice) {
                final ReferenceLocation location = decodeAbsoluteLocation(absoluteLocation);
                location.compressedByteAmountSincePreviousLocation=compressedByteAmountSincePreviousLocation;
                result.add(location);
                lastFileOffsetPushed = offsetInEntriesFile;
            }
            i++;
        }
        return result;
    }

    /**
     * Returns a sample of locations covered by this alignment.
     *
     * @param modulo Modulo to avoid sampling every position in the genome.
     * @return A set of positions that do occur in the genome, rounded to the specified modulo value (absoluteLocation-(absoluteLocation % modulo)).
     *         * @throws IOException
     */
    public ObjectList<ReferenceLocation> getLocations(int modulo) throws IOException {
        assert isHeaderLoaded() : "header must be loaded to query locations.";
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
            if (absoluteLocation > offset - 1) {

                absoluteLocation -= offset;
                break;
            }
        }

        return new ReferenceLocation(referenceIndex, (int) absoluteLocation);
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

    @Override
    public long getStartByteOffset(final int startReferenceIndex, final int startPosition) {
        //  System.out.printf("start target: %d position: %d %n",startReferenceIndex, startPosition);
        return getByteOffset(startReferenceIndex, startPosition, 0);
    }

    @Override
    public boolean getQueryIndicesWerePermuted() {
        return queryIndicesWerePermuted;
    }

    @Override
    public boolean hasQueryIndexOccurrences() {
        return hasQueryIndexOccurrences;
    }

    @Override
    public boolean hasAmbiguity() {
        return hasAmbiguity;
    }

    @Override
    public long getEndByteOffset(final int startReferenceIndex, final int startPosition, final int endReferenceIndex, final int endPosition) {
        final long startByteOffset = getByteOffset(startReferenceIndex, startPosition, 0);
        long endByteOffset = startByteOffset;
        int i = 1;
        while (endByteOffset == startByteOffset) {
            endByteOffset = getByteOffset(endReferenceIndex, endPosition + 1, i);
            //    System.out.println("i="+i);
            ++i;
        }
        //  System.out.println(endByteOffset);
        return endByteOffset;
    }

    @Override
    public ReadOriginInfo getReadOriginInfo() {
        return new ReadOriginInfo(readOriginInfoList);
    }

    /**
     * Return the version of Goby that created this alignment. If the alignment did not store a version number explicitely,
     * the string "1.9.5-" is returned to represent all versions of Goby released before Goby 1.9.6.
     *
     * @return A string in the format 1.9.5, or a timestamp YYYYMMDDHHMMSS, where YYYY is the year, MM is the month, DD is the day, HH is the hour, MM the minutes and SS the seconds (time of compilation) when the alignment was created with a development version of Goby.
     */

    public String getGobyVersion() {
        assert isHeaderLoaded() : "header must be loaded to query Goby version.";
        return gobyVersion == null || "".equals(gobyVersion) ? "1.9.5-" : gobyVersion;
    }
}
