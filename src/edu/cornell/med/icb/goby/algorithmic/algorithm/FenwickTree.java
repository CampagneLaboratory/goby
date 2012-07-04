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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import java.io.Serializable;

/**
 * Fenwick tree implementation.
 */
public class FenwickTree implements Serializable {
    private static final long serialVersionUID = -7830133715336861385L;
    private int n;
    private long totalCount;

    public int size() {
        return n;
    }

    public FenwickTree(int n) {
        this.cumCount = new long[n + 1];
        this.n = n;
    }

    /* The following methods implement a Fenwick tree. See http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=binaryIndexedTrees#reada
    * */
    long[] cumCount;

    /**
     * Increment the count of an element.
     *
     * @param index index of the element.
     */
    public void incrementCount(int index) {
        ++index;
        ++totalCount;
        while (index <= n) {
            cumCount[index]++;
            index += index & -index; // By chance, this gives the right next index 8^).
        }
    }

    /**
     * Get the cumulative count for elements between [0-index].
     *
     * @param index index of the element.
     * @return count for element at index.
     */
    public long getCumulativeCount(int index) {
        if (index >= cumCount.length - 1) {
            // past the capacity of the array is all the counts we have seen:
            return totalCount;
        }
        int count = 0;
        index++;
        while (index != 0) {
            count += cumCount[index];
            index = index & index - 1; // This cancels out the least nonzero bit.
        }

        return count;
    }

    /**
     * Get the cumulative count over all the elements. This is exactly the number of times increment count has
     * been called.
     *
     * @return
     */
    public long getTotalCount() {
        return totalCount;
    }


    /**
     * Retrieve the symbol index that has cumulative count equal to the countQueried argument.
     *
     * @param cumulativeCountQueried count for symbol that we are looking for.
     * @return the countQueried, or -1 if the count does not appear in the tree.
     */
    public int find(final int cumulativeCountQueried, final answer result) {
        int start = 0;
        int end = n;
        int middle = start + end >> 1;
        result.cumulativeCount = -1;

        while (middle >= 0 && middle < end) {
            if (end == start) {
                //     System.out.printf("Returning from findX with countQueried=%d %n", start - 1);
                result.symbolIndex = start - 1;
            }
            //   result.cumulativeCount = getCumulativeCount(middle);
            long count;
            long countIndexPlus1 = 0;
            int index = middle;

            if (middle >= cumCount.length - 1) {
                // past the capacity of the array is all the counts we have seen:
                count = totalCount;
                countIndexPlus1 = totalCount;
            } else {
                count = 0;
                index++;
                int indexPlusZero = index;
                int indexPlus1 = index + 1;
                final int mask = index & indexPlus1;
                int commonIndex = mask;
                while (commonIndex != 0) {
                    final long ccount = cumCount[index];
                    count += ccount;
                    countIndexPlus1 += ccount;
                    commonIndex = commonIndex & commonIndex - 1; // This cancels out the least nonzero bit.
                }
                index = indexPlusZero & ~mask; // the part of indexPlusZero that is not common with indexPlus1
                while (index!=0) {
                    count += cumCount[index];
                    index = index & index - 1; // This cancels out the least nonzero bit.
                }
                index = indexPlus1;
                while (index >indexPlusZero && index <cumCount.length) {
                    countIndexPlus1 += cumCount[index];
                    index = index & index - 1; // This cancels out the least nonzero bit.
                }
            }
            result.cumulativeCount = count;
            result.nextLargerCumulativeCount = countIndexPlus1;
            // System.out.printf("start=%d middle=%d end=%d count=%d %n", start, middle, end, count);
            if (cumulativeCountQueried < result.cumulativeCount) {
                end = middle;
                middle = start + end >> 1;
            } else if (cumulativeCountQueried > result.cumulativeCount) {
                start = middle + 1;
                middle = start + end >> 1;

            } else {
                // lastCount = count;
                // lastIndex = middle;
                result.symbolIndex = middle;
                break;
            }
        }
        if (result.cumulativeCount == cumulativeCountQueried) {
            return result.symbolIndex;
        } else {
            // the queried count does not exist in this tree.
            return -1;
        }

    }

    public answer createAnswer() {
        return new answer();
    }

    public class answer {
        int symbolIndex = -1;
        long cumulativeCount = -1;
        long nextLargerCumulativeCount = -1;
    }
}
