/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
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

package edu.cornell.med.icb.goby.alignments.perms;

import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import it.unimi.dsi.fastutil.ints.*;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.util.BitSet;

/**
 * Builds a permutation of query indices to small values.
 *
 * @author Fabien Campagne
 *         Date: 3/5/12
 *         Time: 5:10 PM
 */
public class QueryIndexPermutation implements QueryIndexPermutationInterface {
    /**
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(QueryIndexPermutation.class);
    @RegisterThis
    public static DynamicOptionClient doc = new DynamicOptionClient(QueryIndexPermutation.class,
            "safe-mode:boolean, when true keeps query indices in memory even when the link appears" +
                    " to point backwards. This can help process some incorrect BAM files where pair-links" +
                    " incorrectly map the mate on the same reference, when it appears on a different chromosome  with " +
                    " a position earlier than the primary read. Please note that this option can consume large amounts of " +
                    " memory and should be used only for problematic BAM input files:false"
    );
    private int smallestIndex = Integer.MAX_VALUE;
    private int biggestSmallIndex = Integer.MIN_VALUE;
    private PermutationWriter permutationWriter;
    private final String basename;
    private int globalQueryMaxOccurences = 1;
    private final Int2IntMap offlinePermutation = new Int2IntLinkedOpenHashMap();
    private static final int MAX_OFFLINE_CAPACITY = 100000;
    private boolean isSafeMode;

    public static DynamicOptionClient doc() {
        return doc;
    }

    public void reset() {
        smallestIndex = Integer.MAX_VALUE;
        biggestSmallIndex = Integer.MIN_VALUE;
        smallIndexCounter = 0;
        queryIndexPermutation.clear();
        queryIndexPermutation.defaultReturnValue(-1);
        timesRequested.defaultReturnValue((byte) 1);
        if (permutationWriter != null) {
            permutationWriter.close();
        }
        permutationWriter = new PermutationWriter(basename);
        isSafeMode = doc().getBoolean("safe-mode");
    }

    public QueryIndexPermutation(String filename) {
        this.basename = AlignmentReaderImpl.getBasename(filename);
        reset();

    }

    @Override
    public Alignments.AlignmentEntry makeSmallIndices(final Alignments.AlignmentEntry entry) {
        final Alignments.AlignmentEntry.Builder merged = Alignments.AlignmentEntry.newBuilder(entry);
        makeSmallIndices(merged);
        return merged.build();
    }

    @Override
    public void makeSmallIndices(final Alignments.AlignmentEntry.Builder entry) {
        final int queryIndex = entry.getQueryIndex();
        /* if (queryIndex == 127177) {
            System.out.println("STOP");
            System.out.flush();
        }*/

        final int maxOccurence = calculateQueryIndexOccurrences(entry);
        /*if (queryIndex == 127177) {
            System.out.println("queryIndex=413876 qio=" + maxOccurence);
        } */

        final int smallIndex = getSmallIndex(queryIndex, maxOccurence);
        entry.setQueryIndex(smallIndex);
        smallestIndex = Math.min(smallestIndex, smallIndex);
        biggestSmallIndex = Math.max(biggestSmallIndex, smallIndex);
        assert smallIndex != -1 : "Query small index must never be negative. Observed for original queryIndex=" +
                queryIndex + " " + entry.build().toString();
    }

    private int calculateQueryIndexOccurrences(final Alignments.AlignmentEntry.Builder entry) {
        if (entry.hasQueryIndexOccurrences()) {

            // if we already know the number of times query index occurs in the genome, use that.
            return entry.getQueryIndexOccurrences();
        } else {
            // all entries have at least one occurrence across the genome (this is why they are in the entries file):
            final int queryIndex = entry.getQueryIndex();
            int queryIndexOccurrences = timesRequested.get(queryIndex) + 1;
            // entries with a paired entry in the future get a +1
            queryIndexOccurrences += entry.hasPairAlignmentLink() && (isForward(entry, entry.getPairAlignmentLink()) || isSafeMode) ? 1 : 0;
            // entries with a spliced link forward get +1
            queryIndexOccurrences += entry.hasSplicedForwardAlignmentLink() ? 1 : 0;
            if (entry.hasPairAlignmentLink()) {
                queryIndexOccurrences = Math.max(queryIndexOccurrences, entry.getPairAlignmentLink().getFragmentIndex() + 1);
            }
            if (entry.hasSplicedForwardAlignmentLink()) {
                queryIndexOccurrences = Math.max(queryIndexOccurrences, entry.getSplicedForwardAlignmentLink().getFragmentIndex() + 1);
            }
            return queryIndexOccurrences;
        }
    }

    // determine if the pair link points forward in genomic orientation. Only needs to deal with sorted alignments
    // since permutations are not written for unsorted.
    private boolean isForward(final Alignments.AlignmentEntry.Builder entry,
                              final Alignments.RelatedAlignmentEntry pairAlignmentLink) {

        final int targetIndex = entry.getTargetIndex();
        final int position = entry.getPosition();
        final int linkedTargetIndex = pairAlignmentLink.getTargetIndex();
        final int linkedPosition = pairAlignmentLink.getPosition();

        if (linkedTargetIndex == targetIndex) {
            return linkedPosition >= position;
        } else {

            return linkedTargetIndex > targetIndex;
        }
    }

    @Override
    public int getSmallestIndex() {
        if (smallestIndex == Integer.MAX_VALUE) {
            return 0;
        }
        return smallestIndex;
    }

    @Override
    public int getBiggestSmallIndex() {
        return biggestSmallIndex;
    }

    @Override
    public final void setSmallestIndex(final int value) {
        smallestIndex = value;
    }

    @Override
    public final void setBiggestSmallIndex(final int value) {
        biggestSmallIndex = value;
    }

    @Override
    public int permutate(final int queryIndex) {
        return permutate(queryIndex, globalQueryMaxOccurences);
    }

    @Override
    public int permutate(final int queryIndex, final int maxQueryIndexOccurrence) {
        final int smallIndex = getSmallIndex(queryIndex, maxQueryIndexOccurrence);
        smallestIndex = Math.min(smallestIndex, smallIndex);
        biggestSmallIndex = Math.max(biggestSmallIndex, smallIndex);
        return smallIndex;
    }

    /**
     * Call this method when you know the queryIndex is currently in the queryIndexPermutation map.
     *
     * @param queryIndex
     * @param maxQueryIndexOccurrence
     * @return
     */
    public int internalDoPerm(final int queryIndex, final int maxQueryIndexOccurrence) {
        final int smallIndex = queryIndexPermutation.get(queryIndex);
        final byte timesSeen = (byte) (timesRequested.get(queryIndex) + 1);
        // decide if we have reached max observations for this query index:
        if (timesSeen >= maxQueryIndexOccurrence) {
            // if yes, remove the index from the map, it will not be asked again.
            queryIndexPermutation.remove(queryIndex);

            pushToPreStorage(queryIndex, smallIndex);
        } else {
            // if not, keep it in the map until requested that many times.
            timesRequested.put(queryIndex, timesSeen);
        }
        queryIndicesAlreadySeen.set(queryIndex);
        return smallIndex;

    }

    private void pushToPreStorage(int queryIndex, int smallIndex) {
        //      System.out.printf("pushing to pre-storage queryIndex=%d smallIndex=%d %n",queryIndex, smallIndex);

        moveIndexToPreOffline(queryIndex, smallIndex);
        if (offlinePermutation.size() > MAX_OFFLINE_CAPACITY) {
            save();
        }
    }

    private int smallIndexCounter = 0;
    private final Int2IntMap queryIndexPermutation = new Int2IntOpenHashMap();
    private final BitSet queryIndicesAlreadySeen = new BitSet();
    private final Int2ByteMap timesRequested = new Int2ByteOpenHashMap();

    private int getSmallIndex(final int queryIndex, final int maxObservations) {
        if (!queryIndicesAlreadySeen.get(queryIndex)) {

            // not seen before, let's associate the next small index for this new query index
            /* final int result = queryIndexPermutation.get(queryIndex);
           assert result==-1 :" the result cannot be different from -1";
           if (result == -1) {
            */
            final int smallIndex = smallIndexCounter++;
            queryIndicesAlreadySeen.set(queryIndex, true);
            if (maxObservations > 1) {

                queryIndexPermutation.put(queryIndex, smallIndex);
                timesRequested.put(queryIndex, (byte) 1);
            } else {
                // if maxObs<=1 we don't need to remember the queryIndex in memory
                pushToPreStorage(queryIndex, smallIndex);
            }

            return smallIndex;
            /* } else {
                return result;
            }*/
        } else {
            // the query index was seen before, and we need to return the small index previously associated with
            // that large index.

            final int smallIndex = internalDoPerm(queryIndex, maxObservations);
            return smallIndex;
        }
        //    return fetchExternal(queryIndex);
    }

    /**
     * Take a query index and associated small index and move to pre-offline (immediate state before write).
     *
     * @param queryIndex
     * @param smallIndex
     */
    private void moveIndexToPreOffline(final int queryIndex, final int smallIndex) {
        offlinePermutation.put(queryIndex, smallIndex);

    }

    public void setPruneLimit(byte limit) {
        globalQueryMaxOccurences = limit;
    }

    @Override
    public void close() {
        // move everything left to pre-offline state:
        for (int queryIndex : queryIndexPermutation.keySet()) {

            moveIndexToPreOffline(queryIndex, queryIndexPermutation.get(queryIndex));
        }

        queryIndexPermutation.clear();
        // now save it:
        save();
        permutationWriter.close();
    }

    /**
     * Return true if the query index is kept in memory with its small index, false otherwise.
     *
     * @param queryIndex
     * @return
     */
    public boolean isInMap(int queryIndex) {
        return queryIndexPermutation.containsKey(queryIndex);
    }

    private void save() {
        if (LOG.isTraceEnabled()) {
            LOG.trace("Saving new permutation chunk ");
        }
        try {
            permutationWriter.append(offlinePermutation);
            offlinePermutation.clear();

        } catch (IOException e) {
            throw new RuntimeException("Unable to write permutation component.", e);
        }


    }

    /**
     * Indicates that the query index is now on disk.
     *
     * @param queryIndex
     * @return
     */
    public boolean isOnDisk(int queryIndex) {
        return !isInMap(queryIndex) && queryIndicesAlreadySeen.get(queryIndex);
    }
}

