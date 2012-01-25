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

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.util.*;

/**
 * Read over a set of alignments. This aligner concatenates entries from the input alignment.
 * Reference sequences must match exactly across the input alignments.
 * Query are assumed to be entirely distinct and will be treated as independent observations (e.g.,
 * reads from multiple independent samples). To this effect, alignment entries read from
 * different input basenames, which would otherwise share an identical query index,
 * are renumbered with distinct query indices.
 *
 * @author Fabien Campagne
 *         Date: May 20, 2009
 *         Time: 5:06:01 PM
 */
public class ConcatAlignmentReader extends AbstractConcatAlignmentReader {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(ConcatAlignmentReader.class);

    protected final AlignmentReader[] readers;
    protected final IntSet readersWithMoreEntries;

    /**
     * One element per reader.
     */
    private final int[] numQueriesPerReader;

    /**
     * One element per reader.
     */
    private final int[] queryIndexOffset;

    protected int activeIndex;
    protected boolean adjustQueryIndices = true;
    private int numberOfAlignedReads;

    /**
     * Construct an alignment reader over a set of alignments.
     * Please note that the constructor access the header of each individual alignment to
     * check reference sequence identity and obtain the number of queries in each input alignment.
     * This version uses adjustQueryIndices as the default true.
     *
     * @param basenames Basenames of the individual alignemnts to combine.
     * @throws IOException If an error occurs reading the header of the alignments.
     */
    public ConcatAlignmentReader(final String... basenames) throws IOException {
        this(new DefaultAlignmentReaderFactory(), true, basenames);
    }

    /**
     * Construct an alignment reader over a set of alignments.
     * Please note that the constructor access the header of each individual alignment to
     * check reference sequence identity and obtain the number of queries in each input alignment.
     * This version uses adjustQueryIndices as the default true.
     *
     * @param basenames          Basenames of the individual alignemnts to combine.
     * @param adjustQueryIndices if we need to adjustQueryIndices
     * @throws IOException If an error occurs reading the header of the alignments.
     */
    public ConcatAlignmentReader(boolean adjustQueryIndices, final String... basenames) throws IOException {
        this(new DefaultAlignmentReaderFactory(), adjustQueryIndices, basenames);
    }

    /**
     * Construct an alignment reader over a set of alignments.
     * Please note that the constructor access the header of each individual alignment to
     * check reference sequence identity and obtain the number of queries in each input alignment.
     *
     * @param alignmentReaderFactory Factory to create new alignmentReaders.
     * @param adjustQueryIndices     if we need to adjustQueryIndices
     * @param basenames              Basenames of the individual alignemnts to combine.
     * @throws IOException If an error occurs reading the header of the alignments.
     */
    public ConcatAlignmentReader(final AlignmentReaderFactory alignmentReaderFactory,
                                 final boolean adjustQueryIndices, final String... basenames) throws IOException {
        super(true, null);
        this.adjustQueryIndices = adjustQueryIndices;
        readers = alignmentReaderFactory.createReaderArray(basenames.length);

        readersWithMoreEntries = new IntArraySet();

        for (int readerIndex = 0; readerIndex < basenames.length; readerIndex++) {
            readers[readerIndex] = alignmentReaderFactory.createReader(basenames[readerIndex]);
            readersWithMoreEntries.add(readerIndex);
            sampleBasenames.add(basenames[readerIndex]);
        }
        numQueriesPerReader = new int[basenames.length];
        queryIndexOffset = new int[basenames.length];
        readHeader();
    }


    /**
     * Construct an alignment reader over a set of alignments.
     * Please note that the constructor access the header of each individual alignment to
     * check reference sequence identity and obtain the number of queries in each input alignment.
     *
     * @param alignmentReaderFactory Factory to create new alignmentReaders.
     * @param adjustQueryIndices     if we need to adjustQueryIndices
     * @param startReferenceIndex    Index of the reference for the start position.
     * @param startPosition          Position on the reference for the start position.
     * @param endReferenceIndex      Index of the reference for the end position.
     * @param endPosition            Position on the reference for the end position.
     * @param basenames              Basenames of the individual alignemnts to combine.
     * @throws IOException If an error occurs reading the header of the alignments.
     */
    public ConcatAlignmentReader(final AlignmentReaderFactory alignmentReaderFactory,
                                 final boolean adjustQueryIndices,
                                 final int startReferenceIndex,
                                 final int startPosition,
                                 final int endReferenceIndex,
                                 final int endPosition,
                                 final String... basenames) throws IOException {
        super(true, null);
        this.adjustQueryIndices = adjustQueryIndices;
        readers = alignmentReaderFactory.createReaderArray(basenames.length);
        readersWithMoreEntries = new IntArraySet();
        int readerIndex = 0;
        for (final String basename : basenames) {
            readers[readerIndex] = alignmentReaderFactory.createReader(basename,
                    startReferenceIndex, startPosition,
                    endReferenceIndex, endPosition);
            readersWithMoreEntries.add(readerIndex);
            sampleBasenames.add(basename);
            readerIndex++;
        }
        numQueriesPerReader = new int[basenames.length];
        queryIndexOffset = new int[basenames.length];
        readHeader();
    }

    /**
     * Read the header of this alignment.
     *
     * @throws java.io.IOException If an error occurs.
     */
    @Override
    public final void readHeader() throws IOException {
        if (!isHeaderLoaded()) {
            final IntSet targetNumbers = new IntArraySet();
            int readerIndex = 0;
            ObjectList<String> alignerNames = new ObjectArrayList<String>();
            ObjectList<String> alignerVersions = new ObjectArrayList<String>();

            numberOfQueries = 0;
            smallestQueryIndex = Integer.MAX_VALUE;
            largestQueryIndex = adjustQueryIndices ? Integer.MIN_VALUE : 0;
            for (final AlignmentReader reader : readers) {
                reader.readHeader();
                String alignerName = reader.getAlignerName();
                String alignerVersion = reader.getAlignerVersion();
                if (!(alignerNames.contains(alignerName) && alignerVersions.contains(alignerVersion))) {
                    alignerNames.add(alignerName);
                    alignerVersions.add(alignerVersion);
                }

                smallestQueryIndex = Math.min(reader.getSmallestSplitQueryIndex(), smallestQueryIndex);
                largestQueryIndex = adjustQueryIndices ?
                        Math.max(largestQueryIndex, 0) + 1 + reader.getLargestSplitQueryIndex() :
                        Math.max(reader.getLargestSplitQueryIndex(), largestQueryIndex);

                targetNumbers.add(reader.getNumberOfTargets());
                final int numQueriesForReader = reader.getNumberOfQueries();
                numQueriesPerReader[readerIndex] = numQueriesForReader;
                numberOfQueries += numQueriesForReader;
                numberOfAlignedReads += reader.getNumberOfAlignedReads();
                readerIndex++;
            }
            alignerName = alignerNames.toString();
            alignerVersion = alignerVersions.toString();
            if (targetNumbers.size() != 1) {
                throw new IllegalArgumentException("The number of targets must match exactly across the input basenames. Found " + targetNumbers.toString());
            } else {
                this.numberOfTargets = targetNumbers.iterator().nextInt();
            }
            targetIdentifiers = new IndexedIdentifier();
            // target information may have more or less targets depending on the reader, but indices must match across
            // all readers:
            boolean error = false;
            for (final AlignmentReader reader : readers) {
                IndexedIdentifier targetIds = reader.getTargetIdentifiers();
                for (MutableString key : targetIds.keySet()) {
                    if (!targetIdentifiers.containsKey(key)) {
                        targetIdentifiers.put(key, targetIds.getInt(key));
                    } else {
                        final int globalValue = targetIdentifiers.getInt(key);
                        final int localValue = targetIds.getInt(key);
                        if (globalValue != localValue) {
                            error = true;
                            LOG.error(
                                    String.format("target indices must match across input alignments. Key %s was found with the distinct values global: %d local %d in alignment %s",
                                            key, globalValue, localValue, reader.basename()));
                        }
                    }
                }
            }
            if (error) {
                throw new RuntimeException("target indices must match across input alignments.");
            }
            targetLengths = new int[targetIdentifiers.size()];
            // keep the maximum length across all readers. We do this to retrieve targetLength over alignments merged
            // from pieces that do not have entries for all target.
            for (int targetIndex = 0; targetIndex < targetIdentifiers.size(); targetIndex++) {
                int maxLength = -1;
                for (final AlignmentReader reader : readers) {
                    final int[] readerLengths = reader.getTargetLength();
                    if (readerLengths != null && readerLengths.length > targetIndex) {
                        maxLength = Math.max(readerLengths[targetIndex], maxLength);
                        targetLengths[targetIndex] = maxLength;
                    }
                }

            }
            // calculate offsets needed to adjustQueryIndices
            for (int i = 0; i < queryIndexOffset.length; i++) {

                queryIndexOffset[i] = adjustQueryIndices ? i == 0 ? 0 : readers[i - 1].getLargestSplitQueryIndex() + 1 : 0;

            }

        }


        setHeaderLoaded(true);
    }


    protected int mergedQueryIndex(final int queryIndex) {
        return queryIndexOffset[activeIndex] + queryIndex;
    }

    /**
     * Iterator over alignment entries.
     *
     * @return an iterator over the alignment entries.
     */
    public final Iterator<Alignments.AlignmentEntry> iterator() {
        return this;
    }

    /**
     * Returns true if the input has more entries.
     *
     * @return true if the input has more entries, false otherwise.
     */
    public boolean hasNext() {
        while (!readersWithMoreEntries.isEmpty()) {
            activeIndex = readersWithMoreEntries.iterator().nextInt();
            final AlignmentReader reader = readers[activeIndex];
            final boolean hasNext = reader.hasNext();
            if (!hasNext) {
                readersWithMoreEntries.remove(activeIndex);
            } else {
                return true;
            }

        }
        return false;
    }


    /**
     * @return The list of aligner names with duplicates removed
     */
    @Override
    public String getAlignerName() {
        return super.getAlignerName();
    }

    /**
     * @return The list of aligner versions with duplicates removed
     */
    @Override
    public String getAlignerVersion() {
        return super.getAlignerVersion();
    }


    /**
     * Returns the next alignment entry from the input stream.
     *
     * @return the alignment read entry from the input stream.
     */
    public Alignments.AlignmentEntry next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        } else {
            final Alignments.AlignmentEntry alignmentEntry = readers[activeIndex].next();
            final int queryIndex = alignmentEntry.getQueryIndex();
            final int newQueryIndex = mergedQueryIndex(queryIndex);


            Alignments.AlignmentEntry.Builder builder = alignmentEntry.newBuilderForType().mergeFrom(alignmentEntry);
            if (adjustQueryIndices && newQueryIndex != queryIndex) {

                builder = builder.setQueryIndex(newQueryIndex);
            }
            if (adjustSampleIndices) {
                builder = builder.setSampleIndex(activeIndex);
            }

            return builder.build();
        }
    }

    /**
     * This operation is not supported by this iterator.
     */
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a reader.");
    }

    /**
     * @deprecated
     */

    @Deprecated
    public void setAdjustQueryIndices(final boolean adjustQueryIndices) {
        throw new UnsupportedOperationException("This operation is unsafe. Set flag through the constructor.");
    }

    /**
     * Obtain statistics about this alignment as a Java property instance.
     *
     * @return statistics about this alignment
     */
    public Properties getStatistics() {
        int index = 1;
        final Properties result = new Properties();
        for (final AlignmentReader reader : this.readers) {
            final Properties localProps = reader.getStatistics();
            for (final Map.Entry<Object, Object> localProp : localProps.entrySet()) {
                result.put("part" + index + "." + localProp.getKey().toString(), localProp.getValue());
            }
            index++;
        }
        return result;
    }

    public int getNumberOfAlignedReads() {
        return numberOfAlignedReads;
    }

    /**
     * Close the underlying readers.
     *
     * @throws IOException if an I/O error occurs
     */
    public void close() throws IOException {
        for (final AlignmentReader reader : readers) {
            reader.close();
        }
    }

    public ObjectList<ReferenceLocation> getLocations(int modulo) throws IOException {
        readHeader();
        ObjectSet<ReferenceLocation> result = new ObjectOpenHashSet<ReferenceLocation>();

        for (AlignmentReader reader : this.readers) {
            result.addAll(reader.getLocations(modulo));
        }
        ObjectList<ReferenceLocation> list = new ObjectArrayList<ReferenceLocation>();
        list.addAll(result);
        Collections.sort(list);
        return list;
    }


}
