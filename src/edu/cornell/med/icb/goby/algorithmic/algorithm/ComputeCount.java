/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.algorithmic.data.Read;
import edu.cornell.med.icb.goby.counts.CountsWriter;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntAVLTreeSet;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntSortedSet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.util.Collections;
import java.util.Map;

/**
 * Data structure and algorithm to compute base-level read coverage histogram over a reference sequence.
 */
public class ComputeCount implements ComputeCountInterface {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(ComputeCount.class);

    /**
     * Number of reads that start ON or BEFORE this position.
     */
    protected final Int2IntMap starts;
    /**
     * Number of reads that end BEFORE position.
     */
    protected final Int2IntMap ends;
    protected final IntList startKeys;
    protected final IntList endKeys;
    /**
     * Used to store bases count in memory.
     */
    protected final Int2IntMap countPerBase;
    /**
     * Sorted keys of countPerBase.
     */
    protected final IntList countKeys;


    protected boolean startPopulateInitialized;

    public ComputeCount() {
        super();
        startKeys = new IntArrayList();
        endKeys = new IntArrayList();
        starts = new Int2IntOpenHashMap();
        ends = new Int2IntOpenHashMap();

        countPerBase = new Int2IntOpenHashMap();
        countKeys = new IntArrayList();
    }

    /**
     * This method must be called before calling the populate method. It initializes data structures.
     */
    public void startPopulating() {
        starts.defaultReturnValue(0);
        ends.defaultReturnValue(0);
        countPerBase.defaultReturnValue(0);
        countPerBase.put(0, 0);
        starts.put(0, 0);
        ends.put(0, 0);
        startPopulateInitialized = true;
    }

    /**
     * Convenience method to populate with a set of reads. This method does not scale, use only for testing
     * with small sets of reads.
     *
     * @param reads
     */
    public final void populate(final ObjectList<Read> reads) {
        startPopulating();

        for (final Read read : reads) {
            final int startIndex = read.start;
            final int endIndex = read.end;
            populate(startIndex, endIndex);
        }
    }

    /**
     * Populate the data structures with a read.
     *
     * @param startIndex Start position of the read ON or BEFORE the startIndex on the reference sequence.
     * @param endIndex   End position of the read BEFORE the endIndex on the reference sequence.
     * @return length of the read.
     */
    public int populate(final int startIndex, final int endIndex) {
        assert startPopulateInitialized : "You must call startPopulating before you can populate.";

        int sval = starts.get(startIndex);
        starts.put(startIndex, ++sval);

        int eval = ends.get(endIndex + 1);
        ends.put(endIndex + 1, ++eval);

        return endIndex - startIndex;
    }

    /**
     * Returns the cumulative start count at position. The number of reads that starts at or before the specified position.
     *
     * @param position position along the reference sequence.
     * @return Cumulative start count.
     */
    protected final int getNumberOfReadsWithStartAt(final int position) {
        return starts.get(position);
    }

    /**
     * Returns the cumulative end count at position. The number of reads that end immediately before
     * the specified position (the cumulative end count).
     *
     * @param position position along the reference sequence.
     * @return The cumulative end count.
     */
    protected final int getNumberOfReadsWithEndAt(final int position) {
        return ends.get(position + 1);
    }

    /**
     * Accumulate start and end counts to produce cumulative start and end count.
     * Pre-condition: the data structures starts and ends must have been populated (see method populate).
     * Post-condition: the data structures starts and end now contain the cumulative start and end counts.
     * It is bad design to reuse the same data structure to store different information, but is there a
     * significant performance advantage in this case?
     * // TODO evaluate and if no performance advantage, refactor to separate cumulative starts, ends and clear the input data structures.
     */
    public void accumulate() {
        LOG.debug("accumulating starts");
        startKeys.addAll(starts.keySet());
        Collections.sort(startKeys);
        endKeys.addAll(ends.keySet());
        Collections.sort(endKeys);
        accumulateOneMap(starts, startKeys);
        LOG.debug("accumulating ends");
        accumulateOneMap(ends, endKeys);
    }

    /**
     * Calculate the cumulative start or end counts. This method is either called with starts and startKeys,
     * or with ends and endKeys. It will calculate the cumulative counts.
     *
     * @param map
     * @param mapKeys
     */
    public final void accumulateOneMap(final Int2IntMap map, final IntList mapKeys) {
        int prevValue = 0;
        for (final int key : mapKeys) {
            final int value = map.get(key);
            final int newValue = prevValue + value;
            map.put(key, newValue);
            prevValue = newValue;
        }
    }

    /**
     * Calculate base counts. Stores the result in the countPerBase and countKey maps.
     */
    public void baseCount() {
        final IntSortedSet joints = new IntAVLTreeSet();
        joints.addAll(starts.keySet());
        joints.addAll(ends.keySet());
        LOG.debug("joints  " + joints);
        final int[] jointsArray = joints.toArray(new int[joints.size()]);
        int prevCount = 0;
        int count;
        LOG.debug("counting");
        int startValue = starts.get(0);
        int endValue = ends.get(0);
        for (int i = 1; i < jointsArray.length; i++) {
            final int curKey = jointsArray[i];
            if (starts.containsKey(curKey)) {
                startValue = starts.get(curKey);
            }
            if (ends.containsKey(curKey)) {
                endValue = ends.get(curKey);
            }
            count = startValue - endValue;
            if (count != prevCount) {
                countPerBase.put(curKey, count);
                prevCount = count;
            }
        }
        countKeys.addAll(countPerBase.keySet());
        Collections.sort(countKeys);
    }

    public Map getCountPerBase() {
        return countPerBase;
    }

    public IntList getCountKeys() {
        return countKeys;
    }

    /**
     * Calculate base counts and write the result to the specified CountsWriter.
     */
    public void baseCount(final CountsWriter writer) throws IOException {
        final IntSortedSet joints = new IntAVLTreeSet();
        joints.addAll(starts.keySet());
        joints.addAll(ends.keySet());

        final int[] jointsArray = joints.toArray(new int[joints.size()]);
        int prevCount = 0;
        int lengthConstant = 0;
        int count;
        int line = 0;
        LOG.debug("counting");
        int startValue = starts.get(0);
        int endValue = ends.get(0);
        for (int i = 1; i < jointsArray.length; i++) {
            if (line % 1000000 == 0) {
                LOG.debug("line " + line);
            }
            line++;
            final int curKey = jointsArray[i];
            final int prevKey = jointsArray[i - 1];
            if (starts.containsKey(curKey)) {
                startValue = starts.get(curKey);
            }
            if (ends.containsKey(curKey)) {
                endValue = ends.get(curKey);
            }
            count = startValue - endValue;
            lengthConstant += curKey - prevKey;
            if (count != prevCount) {
                writer.appendCount(prevCount, lengthConstant);
                prevCount = count;
                lengthConstant = 0;
            }
        }
        writer.close();
    }


    /**
     * Returns the total number of counts on the reference sequence.
     * @return the number of counts on the reference
     */
    public int totalCountOnReference() {
        return starts.get(startKeys.getInt(startKeys.size() - 1));
    }

    /**
     * This implementation ignores strand, but some sub-classes need this information.
     *
     * @param startPosition Start position of a read.
     * @param endPosition   End position of a read.
     * @param forwardStrand True when the read matches the forward strand, false otherwise.
     */

    public void populate(final int startPosition, final int endPosition, final boolean forwardStrand) {
        populate(startPosition, endPosition);
    }

      /**
     * This implementation ignores strand and queryIndex, but some sub-classes need this information.
     *
     * @param startPosition Start position of a read.
     * @param endPosition   End position of a read.
     * @param forwardStrand True when the read matches the forward strand, false otherwise.
     */

    public void populate(final int startPosition, final int endPosition, final boolean forwardStrand, final int queryIndex) {
        populate(startPosition, endPosition);
    }
}
