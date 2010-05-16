/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import edu.cornell.med.icb.goby.algorithmic.data.ReadWithIndex;
import edu.cornell.med.icb.goby.counts.CountsWriter;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.util.Collections;

/**
 * Data structure and algorithm to compute base-level read coverage histogram over a reference sequence.
 */
public class ComputeWeightCount implements AbstractCount {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(ComputeWeightCount.class);

    /**
     * Cumulative weights of reads that start ON or BEFORE this position.
     */
    protected final Int2DoubleMap starts;
    /**
     * Cumulative weights of reads that end BEFORE position.
     */
    protected final Int2DoubleMap ends;
    protected final IntList startKeys;
    protected final IntList endKeys;
    /**
     * Used to store the weight of each base in memory.
     */
    protected final Int2FloatMap countPerBase;
    /**
     * Sorted keys of countPerBase.
     */
    protected final IntList countKeys;

    protected int fixedLength;
    protected boolean isFixedLength;
    protected boolean startPopulateInitialized;
    private FloatArrayList weights;

    public ComputeWeightCount(FloatArrayList weights) {
        super();
        startKeys = new IntArrayList();
        endKeys = new IntArrayList();
        starts = new Int2DoubleOpenHashMap();
        ends = new Int2DoubleOpenHashMap();
        fixedLength = -1;   // not fixed Length
        countPerBase = new Int2FloatOpenHashMap();
        countKeys = new IntArrayList();
        this.weights = weights;
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
    public final void populate(final ObjectList<ReadWithIndex> reads) {
        startPopulating();

        fixedLength = reads.get(0).end - reads.get(0).start; //first read length
        for (final ReadWithIndex read : reads) {
            final int startIndex = read.start;
            final int endIndex = read.end;
            populate(startIndex, endIndex, read.readIndex);
        }
    }

    /**
     * Populate the data structures with a read.
     *
     * @param startIndex Start position of the read ON or BEFORE the startIndex on the reference sequence.
     * @param endIndex   End position of the read BEFORE the endIndex on the reference sequence.
     * @param readIndex  index of the read for the alignment entry being processed.
     * @return length of the read.
     */
    public double populate(final int startIndex, final int endIndex, int readIndex) {
        assert startPopulateInitialized : "You must call startPopulating before you can populate.";

        double sval = starts.get(startIndex);
        sval += weights.get(readIndex);
        starts.put(startIndex, sval);

        double eval = ends.get(endIndex + 1);
        eval += weights.get(readIndex);
        ends.put(endIndex + 1, eval);

        return endIndex - startIndex;
    }

    /**
     * Returns the cumulative start count at position. The number of reads that starts at or before the specified position.
     *
     * @param position position along the reference sequence.
     * @return Cumulative start count.
     */
    protected final double getNumberOfReadsWithStartAt(final int position) {
        return starts.get(position);
    }

    /**
     * Returns the cumulative end count at position. The number of reads that end immediately before
     * the specified position (the cumulative end count).
     *
     * @param position position along the reference sequence.
     * @return The cumulative end count.
     */
    protected final double getNumberOfReadsWithEndAt(final int position) {
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
    public final void accumulateOneMap(final Int2DoubleMap map, final IntList mapKeys) {
        double prevValue = 0;
        for (final int key : mapKeys) {
            final double value = map.get(key);
            final double newValue = prevValue + value;
            map.put(key, newValue);
            prevValue = newValue;
        }
    }

    /**
     * Calculate reweighted base counts. Stores the result in the countPerBase and countKey maps.
     */
    public void baseCount() {
        final IntSortedSet joints = new IntAVLTreeSet();
        joints.addAll(starts.keySet());
        joints.addAll(ends.keySet());
        LOG.debug("joints  " + joints);
        final int[] jointsArray = joints.toArray(new int[joints.size()]);
        double prevCount = 0;
        double count;
        LOG.debug("counting");
        double startValue = starts.get(0);
        double endValue = ends.get(0);
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
                countPerBase.put(curKey, (float)count);
                prevCount = count;
            }
        }
        countKeys.addAll(countPerBase.keySet());
        Collections.sort(countKeys);
    }

    /**
     * Calculate base counts and write the result to the specified CountsWriter.
     */
    public void baseCount(final CountsWriter writer) throws IOException {
        final IntSortedSet joints = new IntAVLTreeSet();
        joints.addAll(starts.keySet());
        joints.addAll(ends.keySet());

        final int[] jointsArray = joints.toArray(new int[joints.size()]);
        double prevCount = 0;
        int lengthConstant = 0;
        double count;
        int line = 0;
        LOG.debug("counting");
        double startValue = starts.get(0);
        double endValue = ends.get(0);
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
                writer.appendCount((int)Math.round(prevCount), lengthConstant);
                prevCount = count;
                lengthConstant = 0;
            }
        }
        writer.close();
    }


    /**
     * Returns the total number of counts on the reference sequence.
     *
     * @return the number of counts on the reference
     */
    public int totalCountOnReference() {
        return (int)Math.round(starts.get(startKeys.getInt(startKeys.size() - 1)));
    }

    /**
     * This implementation ignores strand, but some sub-classes need this information.
     *
     * @param startPosition Start position of a read.
     * @param endPosition   End position of a read.
     * @param forwardStrand True when the read matches the forward strand, false otherwise.
     */

    public void populate(final int startPosition, final int endPosition, final boolean forwardStrand) {
       throw new UnsupportedOperationException("This implementation does not support strand specific populate.");
    }

    public Int2FloatMap getCountPerBase() {
        return countPerBase;
    }

    public IntList getCountKeys() {

        return countKeys;
    }
}