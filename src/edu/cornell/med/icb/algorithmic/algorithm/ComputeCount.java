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

package edu.cornell.med.icb.algorithmic.algorithm;

import edu.cornell.med.icb.algorithmic.data.Read;
import edu.cornell.med.icb.counts.CountsWriter;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntAVLTreeSet;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntSortedSet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Map;

public class ComputeCount {
    /**
     * Number of reads that start ON or BEFORE this position.
     */
    protected Int2IntMap starts;
    /**
     * Number of reads that end BEFORE position.
     */
    protected Int2IntMap ends;
    protected IntList startKeys;
    protected IntList endKeys;
    /**
     * Used to store bases count in memory.
     */
    protected Int2IntMap countPerBase;
    /**
     * Sorted keys of countPerBase.
     */
    protected IntList countKeys;

    protected int fixedLength;
    protected boolean isFixedLength;
    private static final Logger LOG = Logger.getLogger(ComputeCount.class);
    protected boolean startPopulateInitialized;

    /**
     * Positive strand starts position.
     */
    public Int2IntMap posStarts;
    /**
     * Positive strand starts position.
     */
    public Int2IntMap negStarts;

    public ComputeCount() {
        super();
        startKeys = new IntArrayList();
        endKeys = new IntArrayList();
        starts = new Int2IntOpenHashMap(); //new Int2IntAVLTreeMap();
        ends = new Int2IntOpenHashMap(); //Int2IntAVLTreeMap();
        fixedLength = -1; //not fixed Length
        countPerBase = new Int2IntOpenHashMap();
        countKeys = new IntArrayList();

        posStarts = new Int2IntOpenHashMap();
        negStarts = new Int2IntOpenHashMap();
        posStarts.defaultReturnValue(0);
        negStarts.defaultReturnValue(0);
    }

    public void printPosStarts(final String filename, final String chroName) throws IOException {
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(filename, true));
            writer.write(">" + chroName + "\n");
            for (final Map.Entry<Integer, Integer> entry : posStarts.entrySet()) {
                writer.write(entry.getKey() + "\t" + entry.getValue() + "\n");
            }
        } finally {
            IOUtils.closeQuietly(writer);
        }
    }

    public void populatePos(final int startIndex) {
        int sval = posStarts.get(startIndex);
        posStarts.put(startIndex, ++sval);
    }

    public void populateNeg(final int startIndex) {
        int sval = negStarts.get(startIndex);
        negStarts.put(startIndex, ++sval);
    }

    public void startPopulating() {
        starts.defaultReturnValue(0);
        ends.defaultReturnValue(0);
        countPerBase.defaultReturnValue(0);
        countPerBase.put(0, 0);
        starts.put(0, 0);
        ends.put(0, 0);
        startPopulateInitialized = true;
    }

    public final void populate(final ObjectList<Read> reads) {
        startPopulating();

        fixedLength = reads.get(0).end - reads.get(0).start; //first read length
        for (final Read read : reads) {
            final int startIndex = read.start;
            final int endIndex = read.end;
            final int readLength = populate(startIndex, endIndex);
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
     * Eeturns the number of reads starts ON or BEFORE the position.
     *
     * @param position
     * @return
     */
    protected final int getNumberOfReadsWithStartAt(final int position) {
        return starts.get(position);
    }

    /**
     * Get the number of reads ends before the position.
     *
     * @param position
     * @return
     */
    protected final int getNumberOfReadsWithEndAt(final int position) {
        return ends.get(position + 1);
    }

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

    public final void accumulateOneMap(final Int2IntMap map, final IntList mapKeys) {
        int prevValue = 0;
        for (final int key : mapKeys) {
            final int value = map.get(key);
            final int newValue = prevValue + value;
            map.put(key, newValue);
            prevValue = newValue;
        }
    }

    public void baseCount() {
        final IntSortedSet joints = new IntAVLTreeSet();
        joints.addAll(starts.keySet());
        joints.addAll(ends.keySet());
        LOG.debug("joints  " + joints);
//        System.out.println("joints  " + joints);
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
                //               System.out.println("Key  "+curKey+" Count " + count);
                countPerBase.put(curKey, count);
                prevCount = count;
            }
        }
//        countPerBase.put(joints.last()+2,0);
//        System.out.println("countPerbase    "+countPerBase);
        countKeys.addAll(countPerBase.keySet());
        Collections.sort(countKeys);
    }

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
     *
     * @return
     */
    public int totalCountOnReference() {
        return starts.get(startKeys.getInt(startKeys.size() - 1));
    }


    /**
     * Calculate base count by reference for the input. Returns position -> count map.
     *
     * @param reads
     * @return result
     */
    public Int2IntMap baseRunByReference(final ObjectList<Read> reads) {
        final Int2IntMap result = new Int2IntOpenHashMap();
        populate(reads);
        accumulate();

        result.put(0, totalCountOnReference());
        return result;
    }

    private int ONE_MEGABYTE = 1000000;

    private String getHeapSize() {
        final long totalMemory = Runtime.getRuntime().totalMemory();
        final long freeMemory = Runtime.getRuntime().freeMemory();
        final long memoryUseMB = (totalMemory - freeMemory) / ONE_MEGABYTE;
        final long memoryAvailMB = totalMemory / ONE_MEGABYTE;
        return String.format("%,d / %,d MB", memoryUseMB, memoryAvailMB);
    }

    /**
     * Calculate base count by reference for the input. Returns position -> count map.
     *
     * @param reads
     * @return result
     */
    public void baseRun(final ObjectList<Read> reads) throws IOException {
        populate(reads);
        if (LOG.isDebugEnabled()) {
            LOG.debug("populate " + getHeapSize());
        }
        accumulate();
        if (LOG.isDebugEnabled()) {
            LOG.debug("accumulate " + getHeapSize());
        }
        baseCount();
    }

    public IntList getReadStartList(final ObjectList<Read> reads) {
        final IntList readStart = new IntArrayList();
        for (int i = 0; i < reads.size(); i++) {
            readStart.add(reads.get(i).start);
        }
        return readStart;
    }

    /**
     * This implementation ignores strand, but some sub-classes need this information.
     *
     * @param startPosition Start position of a read.
     * @param endPosition   End position of a read.
     * @param forwardStrand  True when the read matches the forward strand, false otherwise.
     */

    public void populate(int startPosition, int endPosition, boolean forwardStrand) {

        populate(startPosition, endPosition);
    }
}
