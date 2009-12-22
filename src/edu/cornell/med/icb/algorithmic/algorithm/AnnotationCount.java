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

import edu.cornell.med.icb.algorithmic.data.Annotation;
import edu.cornell.med.icb.algorithmic.data.Read;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Collections;

public class AnnotationCount {
    public ObjectList<Read> reads;
    public IntList readStart;
    public ComputeCount baseCounter;

    public AnnotationCount() {
        super();
        baseCounter = new ComputeCount();
        readStart = new IntArrayList();
        reads = new ObjectArrayList<Read>();
    }

    public final void populate(final int startPosition, final int endPosition) {
        final Read read = new Read(startPosition, endPosition);
        reads.add(read);
        readStart.add(startPosition);
        baseCounter.populate(startPosition, endPosition);
    }

    public void sortReads() {
        Collections.sort(reads, new Read.ReadSortByStart());
        Collections.sort(readStart);
    }

    /**
     * For hasmap map with keylist key, find the key for index.
     * if index is in map, return the key, otherwise return the key immediately less than index
     *
     * @param searchKey the key to look for
     * @param keyList the list to search
     * @return the index on the keyList for the searchKey postion or immediate previous one
     */
    public final int getIndex(final int searchKey, final IntList keyList) {
        //  if (map.containsKey(index)) return index;
        final int x = Collections.binarySearch(keyList, searchKey); //fails to find index returns position+1
        int index = x < 0 ? -x - 2 : x; //the true index to
        index = (index < 0) ? 0 : index;
        return index;
    }

    /**
     * Returns the value at index at the map either starts or ends.
     *
     * @param searchKey the position on chromosome to get value
     * @param keyList e.g. startKeys or endKeys
     * @param map e.g. starts or ends
     * @return, the count on position index on chromosome
     */
    public final int getValue(final int searchKey, final IntList keyList, final Int2IntMap map) {
        if (map.containsKey(searchKey)) {
            return map.get(searchKey);
        }
        final int index = getIndex(searchKey, keyList);
        return map.get(keyList.get(index));
    }

    /**
     * Returns the average read coverage per base for a segment, or the sequencing depth.
     * requires baseRun to generate the countPerBase hashmap
     *
     * @param start and end of a segment
     * @return average count per base on this segment
     */
    public final float averageReadsPerPosition(final int start, final int end) {
        long sum = 0;

        final int startIndex = getIndex(start, baseCounter.countKeys);
        final int startKey = baseCounter.countKeys.get(startIndex);
        final int startCount = baseCounter.countPerBase.get(startKey);

        //    System.out.println("start   "+start +"  start index   "+startIndex +"   start Key   "+startKey);
        final long startOverCountArea = (start - startKey) * startCount;

        final int maxIndex = baseCounter.countKeys.size() - 1;
        int endIndex = maxIndex;
        int index = startIndex;
        while (index < maxIndex) {
            final int key = baseCounter.countKeys.get(index);
            final int count = baseCounter.countPerBase.get(key);
            final int nextKey = baseCounter.countKeys.get(index + 1);
            //      System.out.println(key+"    nextkey "+nextKey+" count   "+count);
            if (nextKey > end) {
                endIndex = index;
                break;
            }
            final long recArea = count * (nextKey - key);
            sum += recArea;
            index++;
        }

        final int endKey = baseCounter.countKeys.get(endIndex);
        final int endCount = baseCounter.countPerBase.get(endKey);
        final long endUnderCountArea = (end - endKey + 1) * endCount;
        //    System.out.println("sum "+sum+" startOver   "+startOverCountArea+"  endUnder    "+endUnderCountArea);
        sum = sum - startOverCountArea + endUnderCountArea;

        final int segmentSize = end - start + 1;
        //     System.out.println("new sum "+sum+" segsize "+segmentSize);
        if (end < start) {
            return 0;
        } else {
            final NumberFormat formatter = new DecimalFormat("#.00");
            return Float.parseFloat(formatter.format((float) sum / segmentSize));
        }
    }

    /**
     * Return the number of annotation (gene) overlapped reads.
     *
     * @param annot an annotation
     * @return number of reads covered on the genes except all reads exclusively in introns
     */
    public int geneExpressionCount(final Annotation annot) {
        int sum = readsOverlapSegmentCount(annot.getStart(), annot.getEnd());
        final int numIntrons = annot.segments.size() - 1;
        for (int k = 0; k < numIntrons; k++) {
            sum -= readsInSegmentCount(annot.segments.get(k).end + 1, annot.segments.get(k + 1).start - 1);
        }
        return sum;
    }

    /**
     * Returns the number of reads exclusively in the segment.
     *
     * @param start and end
     * @return the number of reads
     */
    public int readsInSegmentCount(final int start, final int end) {
        final int n = reads.size();
        int i = getIndex(start, readStart);
        int count = 0;

        while (i < n) {
            final Read read = reads.get(i);
            if (read.start >= start && read.end <= end) {
                count++;
            } else if (read.start > end) {
                break;
            }
            i++;
        }
        return count;
    }

    /**
     * Returns the number of reads with any overlap on a given segment.
     *
     * @param start
     * @param end
     * @return the number of reads
     */
    public int readsOverlapSegmentCount(final int start, final int end) {
        return getValue(end, baseCounter.startKeys, baseCounter.starts) -
                getValue(start, baseCounter.endKeys, baseCounter.ends);
    }

    //    public Object2IntMap<String> annotationRun(Object2ObjectMap<String, Annotation> annots, ObjectList<Read> reads) {
//        populate(reads);
//        accumulate();
//        Collections.sort(reads, new Read.ReadSortByStart());
//        //sorted read starts list
//        IntList readStart = getReadStartList(reads);
//        Object2IntMap<String> result = new Object2IntOpenHashMap<String>();
//        for (Object2ObjectMap.Entry<String, Annotation> entry : annots.object2ObjectEntrySet()) {
//            Annotation annot = entry.getValue();
//            int exonCount = exonCounts(annot);
//            int intronCount;
//            if (fixedLength <= 0) {
//                intronCount = intronCounts(annot, reads, readStart);
//            } else {
//                intronCount = intronCountsFixedLength(annot, fixedLength);
//            }
//            result.put(annot.id, exonCount - intronCount);
//        }
//        return result;
//    }
    //require sorted reads and annots
//    public Object2IntMap<String> completeOverlapping(ObjectList<Segment> annots, ObjectList<Read> reads) {
//        PriorityQueue<Read> que = new PriorityQueue<Read>(10, new Read.ReadSortByEnd());
//        int m = annots.size();
//        int n = reads.size();
//        Object2IntMap <String> hash = new Object2IntOpenHashMap<String>();
//        int i = 0;
//        for (int k = 0; k < m; k++) {
//            Segment annot = annots.get(k);
//            while (!que.isEmpty() && que.peek().end < annot.end)
//                que.remove();
//            while (i < n && reads.get(i).start <= annot.start) {
//                if (reads.get(i).end >= annot.end)
//                    que.add(reads.get(i)); //insert in order
//                i++;
//            }
//
//            hash.put(annot.id, que.size());
//        }
//        return hash;
//    }

    //require sorted reads and annots

//
//    public int intronCountsFixedLength(Annotation annot, int fixedLength) {
//        int sum = 0;
//        int numIntrons = annot.segments.size() - 1;
//        for (int k = 0; k < numIntrons; k++) {
//            int intronStart = annot.segments.get(k).end + 1;
//            int intronEnd = annot.segments.get(k + 1).start - 1;
//            int count = getValue(intronEnd - fixedLength, startKeys, starts) - getValue(intronStart, startKeys, starts);
//            sum += count;
//        }
//        return sum;
//    }
//
}
