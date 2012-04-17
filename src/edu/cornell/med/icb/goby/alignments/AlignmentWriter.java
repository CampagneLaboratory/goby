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

import edu.cornell.med.icb.goby.alignments.perms.NoOpPermutation;
import edu.cornell.med.icb.goby.alignments.perms.QueryIndexPermutation;
import edu.cornell.med.icb.goby.alignments.perms.QueryIndexPermutationInterface;
import edu.cornell.med.icb.goby.compression.MessageChunksWriter;
import edu.cornell.med.icb.goby.modes.GobyDriver;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.*;
import java.util.Map;
import java.util.Properties;
import java.util.zip.GZIPOutputStream;

/**
 * This class write alignments in a Protocol Buffer, compressed, binary and splitable format.
 * See Alignements.proto for the specification of this format.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 5:53:50 PM
 */
public class AlignmentWriter implements Closeable {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(AlignmentWriter.class);

    private final Alignments.AlignmentCollection.Builder collectionBuilder;
    private final MessageChunksWriter entriesChunkWriter;
    private IndexedIdentifier queryIdentifiers;
    private IndexedIdentifier targetIdentifiers;
    private boolean headerWritten;
    private final GZIPOutputStream headerOutput;
    private boolean entriesHaveQueryLength;
    private boolean entriesHaveQueryIndexOccurrences = true;
    private boolean allReadQualityScores = true;
    @RegisterThis
    public static DynamicOptionClient doc = new DynamicOptionClient(AlignmentWriter.class,
            "permutate-query-indices:boolean, when true permutates query indices to small values (improves compression, but looses the ability to track alignments back to reads):false"
    );
    private ObjectArrayList<Alignments.ReadOriginInfo.Builder> readOriginInfoBuilderList;
    private boolean entriesHaveAmbiguity = true;

    public static DynamicOptionClient doc() {
        return doc;
    }

    /**
     * Details about aligner.
     */
    private String alignerName;
    private String alignerVersion;
    /**
     * The version of Goby that created this alignment file.
     */
    private String gobyVersion;
    /**
     * Set of query lengths to determine whether or not they are unique.
     */
    private final IntSet uniqueQueryLengths = new IntOpenHashSet();

    /**
     * Length of each target sequence.
     */
    private int[] targetLengths;

    private String[] queryIdentifiersArray;
    private String[] targetIdentifiersArray;
    private int maxTargetIndex = -1;
    private int minQueryIndex = Integer.MAX_VALUE;
    private int maxQueryIndex = Integer.MIN_VALUE;
    private int actualNumberOfQueries = Integer.MIN_VALUE;
    private final Properties stats;
    private boolean statsWritten;
    private final FileWriter statsWriter;
    private final String basename;
    private int numberOfAlignedReads;
    private int constantQueryLength;
    private boolean isConstantQueryLength;
    private Alignments.AlignmentEntry.Builder newEntry;
    private final FileOutputStream alignmentEntries;

    private boolean sortedState;

    // data structures to build index:
    private int previousChunkOffset = -1;
    private int firstTargetIndexInChunk;
    private boolean firstEntryInChunk = true;
    private int firstPositionInChunk;
    private final LongArrayList indexOffsets = new LongArrayList();
    private final LongArrayList indexAbsolutePositions = new LongArrayList();
    private boolean indexWritten;
    private long[] targetPositionOffsets;


    private QueryIndexPermutationInterface permutator;
    private boolean queryIndicesWerePermuted;

    public AlignmentWriter(final String outputBasename) throws IOException {
        alignmentEntries = new FileOutputStream(outputBasename + ".entries");
        headerOutput = new GZIPOutputStream(new FileOutputStream(outputBasename + ".header"));
        statsWriter = new FileWriter(outputBasename + ".stats");
        this.basename = outputBasename;
        collectionBuilder = Alignments.AlignmentCollection.newBuilder();
        entriesChunkWriter = new MessageChunksWriter(alignmentEntries);
        entriesChunkWriter.setParser(new AlignmentCollectionHandler());
        newEntry = Alignments.AlignmentEntry.newBuilder();
        queryIdentifiers = new IndexedIdentifier();
        targetIdentifiers = new IndexedIdentifier();
        stats = new Properties();
        // we assume stats were written until a client puts stats in this writer.
        statsWritten = true;
        setPermutation(doc.getBoolean("permutate-query-indices"));

    }

    /**
     * Indicate whether query indices are small indices. When true, create a PermutationReader to reconstitute
     * original query indices from the stored small indices.
     *
     * @param state True or False.
     */
    public void setPermutation(boolean state) {
        if (state) {
            permutator = new QueryIndexPermutation(basename);
            queryIndicesWerePermuted = true;
        } else {
            permutator = new NoOpPermutation();
            queryIndicesWerePermuted = false;
        }
    }

    public void setSorted(final boolean sortedState) {
        this.sortedState = sortedState;
        if (sortedState) {
            if (targetPositionOffsets == null) {

                final String s = "Indexing sorted alignments requires"
                        + " knowning target lengths before entries are appended. setTargetLength"
                        + " must be called before setSorted(true).";
                LOG.error(s);
                throw new UnsupportedOperationException(s);
            }
        } else {
            // unsorted alignments are never permuted since the original query indices are already monotonically
            // increasing in each chunk.
            setPermutation(false);
        }
    }

    /**
     * Set the query index for the next alignment extry.
     *
     * @param queryIndex The query index of the next alignment entry to append
     * @deprecated use
     *             {@link #appendEntry(edu.cornell.med.icb.goby.alignments.Alignments.AlignmentEntry)}.
     */
    @Deprecated
    public final void setQueryIndex(final int queryIndex) {
        newEntry.setQueryIndex(queryIndex);
    }

    /**
     * Set the target index for the next alignment extry.
     *
     * @param targetIndex The target index of the next alignment entry to append
     * @deprecated use
     *             {@link #appendEntry(edu.cornell.med.icb.goby.alignments.Alignments.AlignmentEntry)}.
     */
    @Deprecated
    public final void setTargetIndex(final int targetIndex) {
        newEntry.setTargetIndex(targetIndex);
    }

    /**
     * Set the target position for the next alignment extry.
     *
     * @param position The target position of the next alignment entry to append
     * @deprecated use
     *             {@link #appendEntry(edu.cornell.med.icb.goby.alignments.Alignments.AlignmentEntry)}.
     */
    @Deprecated
    public final void setTargetPosition(final int position) {
        newEntry.setPosition(position);
    }

    /**
     * Set the score for the next alignment extry.
     *
     * @param score The score of the next alignment entry to append
     * @deprecated use
     *             {@link #appendEntry(edu.cornell.med.icb.goby.alignments.Alignments.AlignmentEntry)}.
     */
    @Deprecated
    public final void setAlignmentScore(final float score) {
        newEntry.setScore(score);
    }

    public void setNumAlignmentEntriesPerChunk(final int numEntriesPerChunk) {
        entriesChunkWriter.setNumEntriesPerChunk(numEntriesPerChunk);
    }

    /**
     * Set fields for the next alignment extry.
     *
     * @param queryIndex           The query index of the next alignment entry to append
     * @param targetIndex          The target index of the next alignment entry to append
     * @param position             The target position of the next alignment entry to append
     * @param score                The score of the next alignment entry to append
     * @param matchesReverseStrand true if the entry matches the reverse strand
     * @deprecated use
     *             {@link #appendEntry(edu.cornell.med.icb.goby.alignments.Alignments.AlignmentEntry)}.
     */
    @Deprecated
    public final void setAlignmentEntry(final int queryIndex, final int targetIndex,
                                        final int position,
                                        final float score,
                                        final boolean matchesReverseStrand, int queryLength) {
        newEntry.setQueryIndex(queryIndex);
        newEntry.setTargetIndex(targetIndex);
        newEntry.setScore(score);
        newEntry.setPosition(position);
        newEntry.setMatchingReverseStrand(matchesReverseStrand);
        newEntry.setMultiplicity(1);
        newEntry.setQueryLength(queryLength);
    }

    /**
     * Obtain the alignment entry that is being prepared. Set values on the entry,
     * then call {@link #appendEntry()}.
     *
     * @return the current alignment entry.
     * @deprecated use
     *             {@link #appendEntry(edu.cornell.med.icb.goby.alignments.Alignments.AlignmentEntry)}.
     */
    @Deprecated
    public Alignments.AlignmentEntry.Builder getAlignmentEntry() {
        return newEntry;
    }

    /**
     * Append the current entry to the file being written.
     *
     * @throws IOException if the entry cannot be written
     * @deprecated use
     *             {@link #appendEntry(edu.cornell.med.icb.goby.alignments.Alignments.AlignmentEntry)}.
     */
    @Deprecated
    public synchronized void appendEntry() throws IOException {
        // update the unique query length set
        uniqueQueryLengths.add(newEntry.getQueryLength());

        maxTargetIndex = Math.max(newEntry.getTargetIndex(), maxTargetIndex);
        permutator.makeSmallIndices(newEntry);
        if (newEntry.hasMultiplicity() && newEntry.getMultiplicity() == 1) {
            newEntry.clearMultiplicity();
        }
        final Alignments.AlignmentEntry builtEntry = newEntry.build();

        this.collectionBuilder.addAlignmentEntries(builtEntry);

        writeIndexEntry(builtEntry);
        newEntry = Alignments.AlignmentEntry.newBuilder();
    }

    private void writeIndexEntry(final Alignments.AlignmentEntryOrBuilder builtEntry) throws IOException {
        // detect when all entries have query-index-occurrences:
        entriesHaveQueryIndexOccurrences &= builtEntry.hasQueryIndexOccurrences();

// detect when all entries have ambiguity:
        entriesHaveAmbiguity &= builtEntry.hasAmbiguity();

        //detect when all entries have read_quality_scores:
        allReadQualityScores &= builtEntry.hasReadQualityScores();

        // detect when one or more entries have query length:
        entriesHaveQueryLength |= builtEntry.hasQueryLength();

        if (firstEntryInChunk) {
            firstTargetIndexInChunk = builtEntry.getTargetIndex();
            firstPositionInChunk = builtEntry.getPosition();
            firstEntryInChunk = false;
        }
        final long currentChunkOffset = entriesChunkWriter.writeAsNeeded(collectionBuilder,
                builtEntry.hasMultiplicity() ? builtEntry.getMultiplicity() : 1);
        // LOG.warn(String.format("#entriesWritten: %d currentChunkOffset: %d previousChunkOffset: %d",
        //        entriesChunkWriter.getTotalEntriesWritten(), currentChunkOffset, previousChunkOffset));
        if (sortedState && entriesChunkWriter.getAppendedInChunk() == 0) {
            // we have just written a new chunk.
            pushIndex(currentChunkOffset, firstTargetIndexInChunk, firstPositionInChunk);
            firstEntryInChunk = true;


        } else {
            firstEntryInChunk = false;

        }

    }

    private void pushIndex(final long startOfChunkOffset, final int firstTargetIndexInChunk, final int firstPositionInChunk) {
        final long newOffset = Math.max(startOfChunkOffset, 0);
        final int size = indexAbsolutePositions.size();
        // remove duplicates because the behavior of binary search is undefined for duplicates:
        /**
         * Keep only the first absolutePosition we encounter and its offset in the file. This is done because if an
         * absolute position repeats at the beginning of several consecutive chunks, we want to keep only the first.
         * Also, binarySearch, which we use to access the indexAbsolutePositions array when the alignment is read
         * has undefined behavior when duplicates exist in the array.
         */
        //      LOG.warn(String.format("INDEX attempting to push offset %d %d %n", firstTargetIndexInChunk, firstPositionInChunk));
        final long codedPosition = recodePosition(firstTargetIndexInChunk, firstPositionInChunk);

        if (size == 0 || codedPosition != indexAbsolutePositions.get(size - 1)) {

            indexOffsets.add(newOffset);
            indexAbsolutePositions.add(codedPosition);
            //    LOG.warn(String.format("INDEX Pushing offset=%d position=%d", newOffset, codedPosition));
        }
    }

    protected long recodePosition(final int firstTargetIndexInChunk, final int firstPositionInChunk) {
        assert firstTargetIndexInChunk < targetPositionOffsets.length : "Target length array must have enough elements to store each possible target index.";
        return targetPositionOffsets[firstTargetIndexInChunk] + firstPositionInChunk;
    }

    /**
     * Append an entry to the file being written. The entry must have been prepared
     * outside this writer.
     *
     * @param entry The entry to append
     * @throws IOException If an error occurs writing this entry.
     */
    public synchronized void appendEntry(Alignments.AlignmentEntry entry) throws IOException {
        if (entry.hasQueryLength()) {
            // update the unique query length set
            uniqueQueryLengths.add(entry.getQueryLength());
        }
        if (entry.hasMultiplicity() && entry.getMultiplicity() == 1) {
            // we remove the multiplicity field, since there is no point in storing the default value:
            newEntry = Alignments.AlignmentEntry.newBuilder(entry);
            newEntry.clearMultiplicity();
            entry = newEntry.build();
        }
        maxTargetIndex = Math.max(entry.getTargetIndex(), maxTargetIndex);
        entry = permutator.makeSmallIndices(entry);
        collectionBuilder.addAlignmentEntries(entry);
        writeIndexEntry(entry);

        numberOfAlignedReads += 1;
    }

    /**
     * Replace the target entry, then append the entry.
     *
     * @param entry
     * @param replacementTargetIndex the replacement targetIndex.
     * @throws IOException
     */
    public void appendEntry(Alignments.AlignmentEntry entry, final int replacementTargetIndex) throws IOException {
        entry = entry.newBuilderForType().mergeFrom(entry).setTargetIndex(replacementTargetIndex).build();
        this.appendEntry(entry);
    }

    /**
     * Set the query length on this entry, then append.
     *
     * @param entry       Entry to write
     * @param queryLength the query length.
     * @throws IOException
     */
    public void appendEntryWithLength(Alignments.AlignmentEntry entry, final int queryLength) throws IOException {
        entry = entry.newBuilderForType().mergeFrom(entry).setQueryLength(queryLength).build();
        this.appendEntry(entry);
    }

    public boolean entriesHaveQueryIndexOccurrences() {
        return entriesHaveQueryIndexOccurrences;
    }

    public boolean isAllReadQualityScores() {

        return allReadQualityScores;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws IOException {

        writeHeader();
        writeStats();

        IOUtils.closeQuietly(headerOutput);
        entriesChunkWriter.close(collectionBuilder);
        if (sortedState && targetPositionOffsets != null) {
            writeIndex();
        }

        IOUtils.closeQuietly(alignmentEntries);
        IOUtils.closeQuietly(statsWriter);
    }

    private void writeIndex() throws IOException {
        if (!indexWritten) {
            // Push the last chunkoffset:
            pushIndex(entriesChunkWriter.getCurrentChunkStartOffset(),
                    firstTargetIndexInChunk, firstPositionInChunk);
            GZIPOutputStream indexOutput = null;
            try {
                indexOutput = new GZIPOutputStream(new FileOutputStream(basename + ".index"));
                final Alignments.AlignmentIndex.Builder indexBuilder = Alignments.AlignmentIndex.newBuilder();
                assert (indexOffsets.size() == indexAbsolutePositions.size()) : "index sizes must be consistent.";
                indexBuilder.addAllOffsets(indexOffsets);
                indexBuilder.addAllAbsolutePositions(indexAbsolutePositions);
                indexBuilder.build().writeTo(indexOutput);
            } finally {
                if (indexOutput != null) indexOutput.close();
                indexWritten = true;
            }
        }
    }

    private synchronized void writeHeader() throws IOException {
        if (!headerWritten) {

            final Alignments.AlignmentHeader.Builder headerBuilder = Alignments.AlignmentHeader.newBuilder();
            // record the version of Goby that created this alignment..
            final String version = VersionUtils.getImplementationVersion(GobyDriver.class);
            headerBuilder.setVersion(version);

            headerBuilder.setLargestSplitQueryIndex(permutator.getBiggestSmallIndex());
            headerBuilder.setSmallestSplitQueryIndex(permutator.getSmallestIndex());
            headerBuilder.setNumberOfTargets(maxTargetIndex + 1);
            headerBuilder.setNumberOfQueries(getNumQueries());
            headerBuilder.setSorted(sortedState);
            headerBuilder.setQueryIndicesWerePermuted(queryIndicesWerePermuted);
            headerBuilder.setQueryIndexOccurrences(entriesHaveQueryIndexOccurrences);
            headerBuilder.setAllReadQualityScores(allReadQualityScores);
            // The Java implementation always indexes an index written in sorted order.
            headerBuilder.setIndexed(sortedState);

            headerBuilder.setQueryNameMapping(getMapping(queryIdentifiers, queryIdentifiersArray));
            headerBuilder.setTargetNameMapping(getMapping(targetIdentifiers, targetIdentifiersArray));
            headerBuilder.setNumberOfAlignedReads(numberOfAlignedReads);
            if (alignerName != null) {
                headerBuilder.setAlignerName(alignerName);
            }
            if (alignerVersion != null) {
                headerBuilder.setAlignerVersion(alignerVersion);
            }
            // determine query lengths are constant (regardless of where they came from)
            if (uniqueQueryLengths.size() == 1) {
                // detected constant read length.
                constantQueryLength = uniqueQueryLengths.iterator().nextInt();
                headerBuilder.setConstantQueryLength(constantQueryLength);
                isConstantQueryLength = true;
            } else {
                constantQueryLength = 0;
                isConstantQueryLength = false;
            }
            if (readOriginInfoBuilderList != null) {
                for (final Alignments.ReadOriginInfo.Builder builder : readOriginInfoBuilderList) {
                    headerBuilder.addReadOrigin(builder);
                }
            }
            headerBuilder.setQueryLengthsStoredInEntries(true);
            headerBuilder.setAmbiguityStoredInEntries(entriesHaveAmbiguity);

            // store target lengths:
            if (targetLengths != null) {
                headerBuilder.addAllTargetLength(IntArrayList.wrap(targetLengths));
            }

            headerBuilder.setVersion(VersionUtils.getImplementationVersion(AlignmentWriter.class));
            headerBuilder.build().writeTo(headerOutput);
            headerWritten = true;
        }
    }

    private synchronized void writeStats() throws IOException {
        if (!statsWritten) {
            stats.put("basename", FilenameUtils.getBaseName(basename));
            stats.put("min.query.index", Integer.toString(permutator.getBiggestSmallIndex()));
            stats.put("max.query.index", Integer.toString(permutator.getSmallestIndex()));
            stats.put("number.of.queries", Integer.toString(getNumQueries()));

            stats.put("basename.full", basename);
            stats.put("number.alignment.entries", Integer.toString(numberOfAlignedReads));
            stats.store(statsWriter, "Statistics for merged alignment. ");

            statsWritten = true;
        }
    }

    /**
     * Provide query identifiers as an array of strings, where queryIndex is the index of the
     * element in the array.
     *
     * @param queryIdentifiersArray Array of query identfiers
     */
    public void setQueryIdentifiersArray(final String[] queryIdentifiersArray) {
        this.queryIdentifiersArray = queryIdentifiersArray;
        maxQueryIndex = queryIdentifiersArray.length - 1;
    }

    /**
     * Provide target identifiers as an array of string, where targetIndex is the index of the
     * element in the array.
     *
     * @param targetIdentifiersArray Array of target identfiers
     */
    public void setTargetIdentifiersArray(final String[] targetIdentifiersArray) {
        this.targetIdentifiersArray = targetIdentifiersArray;
        maxTargetIndex = targetIdentifiersArray.length - 1;
    }

    private Alignments.IdentifierMapping getMapping(final IndexedIdentifier identifiers,
                                                    final String[] targetIdentifiersArray) {
        final Alignments.IdentifierMapping.Builder result = Alignments.IdentifierMapping.newBuilder();

        final ObjectList<Alignments.IdentifierInfo> mappings =
                new ObjectArrayList<Alignments.IdentifierInfo>();

        if (targetIdentifiersArray == null) {
            for (final MutableString id : identifiers.keySet()) {

                mappings.add(Alignments.IdentifierInfo.newBuilder()
                        .setName(id.toString())
                        .setIndex(identifiers.get(id))
                        .build());

            }
        } else {
            // prefer to use the array when available, for speed:

            for (int index = 0; index < targetIdentifiersArray.length; index++) {
                mappings.add(Alignments.IdentifierInfo.newBuilder()
                        .setName(targetIdentifiersArray[index])
                        .setIndex(index)
                        .build());

            }

        }

        result.addAllMappings(mappings);
        return result.build();
    }

    public void printStats(final PrintStream out) {
        this.entriesChunkWriter.printStats(out);
        out.println("Min query index: " + minQueryIndex);
        out.println("Max query index: " + maxQueryIndex);
        out.println("Number of queries: " + getNumQueries());
        out.println("Number of targets: " + (maxTargetIndex + 1));
    }


    public void setQueryIdentifiers(final IndexedIdentifier queryIdentifiers) {
        this.queryIdentifiers = queryIdentifiers;
        for (final int index : queryIdentifiers.values()) {
            permutator.permutate(index, 2);
        }
    }

    public void setTargetIdentifiers(final IndexedIdentifier targetIdentifiers) {
        this.targetIdentifiers = targetIdentifiers;
        for (final int index : targetIdentifiers.values()) {
            maxTargetIndex = Math.max(maxTargetIndex, index);
        }
    }


    public void setTargetLengths(final int[] targetLengths) {
        assert targetLengths != null : "Target lengths cannot be null.";
        assert targetLengths.length > maxTargetIndex
                : "The number of elements of targetLength is too small to accommodate targetIndex="
                + maxTargetIndex;

        // calculate the coding offset for each target index. This information will be used by recode
        targetPositionOffsets = new long[targetLengths.length];
        if (targetLengths.length > 0) {
            targetPositionOffsets[0] = 0;
            for (int targetIndex = 1; targetIndex < targetLengths.length; targetIndex++) {
                targetPositionOffsets[targetIndex] =
                        targetLengths[targetIndex - 1] +
                                targetPositionOffsets[targetIndex - 1];

            }
        }
        this.targetLengths = targetLengths;
    }

    /**
     * Set the total number of queries.
     *
     * @param numQueries The number of query sequences.
     */
    public void setNumQueries(final int numQueries) {
        actualNumberOfQueries = numQueries;
        maxQueryIndex = numQueries - 1;
    }

    /**
     * Get the total number of queries.
     *
     * @return The number of query sequences.
     */
    public int getNumQueries() {
        if (actualNumberOfQueries != Integer.MIN_VALUE) {
            return actualNumberOfQueries;
        } else {

            return permutator.getBiggestSmallIndex() - permutator.getSmallestIndex() + 1;
        }
    }

    /**
     * Set the total number of targets.
     *
     * @param numTargets The number of target sequences.
     */
    public void setNumTargets(final int numTargets) {
        maxTargetIndex = numTargets - 1;
    }

    public void putStatistic(final String description, final String value) {
        statsWritten = false;
        stats.put(description, value);
    }

    public void putStatistic(final String description, final double value) {
        putStatistic(description, String.format("%3.3g", value));
    }

    public void putStatistic(final String description, final int value) {
        putStatistic(description, String.format("%d", value));
    }

    public void setStatistics(final Properties statistics) {
        for (final Map.Entry<Object, Object> statistic : statistics.entrySet()) {
            putStatistic(statistic.getKey().toString(), statistic.getValue().toString());
        }
    }

    public void setSmallestSplitQueryIndex(final int smallestQueryIndex) {
        permutator.setSmallestIndex(smallestQueryIndex);
    }

    public void setLargestSplitQueryIndex(final int largestQueryIndex) {
        permutator.setBiggestSmallIndex(largestQueryIndex);
    }


    public String getAlignerVersion() {
        return alignerVersion;
    }

    public void setAlignerVersion(String alignerVersion) {
        this.alignerVersion = alignerVersion;
    }

    public String getAlignerName() {
        return alignerName;
    }

    public void setAlignerName(String alignerName) {
        this.alignerName = alignerName;
    }


    public void setReadOriginInfo(ObjectArrayList<Alignments.ReadOriginInfo.Builder> readOriginInfoBuilderList) {
        this.readOriginInfoBuilderList = readOriginInfoBuilderList;
    }
}
