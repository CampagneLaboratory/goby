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
    public int getCumulativeCount(int index) {
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
     * @return
     */
    public long getTotalCount() {
        return totalCount;
    }
}
