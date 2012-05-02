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

package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

/**
 * Class that populates a cumulative count array of a predetermined size
 * Facilitates adding new counts to the rightmost element in the array without changing
 *
 * @Author: Nyasha Chambwe
 * Date: 3/1/12
 * Time: 4:14 PM
 */
public class CumulativeArrayKeeperFast {


    /**
     * Cumulative count array
     */
    int[] cumC;
    /*
    *Number of elements that are in the array
    * **/
    int capacityMonitor;
    int indexAtLast;
    /*
    Maximum size of the cumulative count array
     */
    private final int capacity;

    public CumulativeArrayKeeperFast(final int n) {
        if (!(n > 0)) {
            throw new IllegalArgumentException("The declared size of the cumulative count array must be at least 1");
        }
        cumC = new int[n];
        capacityMonitor = 0;
        capacity = n;
        indexAtLast = capacity - 1;
    }

    public void addToRight(final int countAtNewSite) {

        if (checkCumulativeCapacityReached()) {
            shift(countAtNewSite);

        } else {
            fillUpCumulativeArray(countAtNewSite);
            capacityMonitor++;
        }
    }

    private void shift(final int countAtNewSite) {
        final int removedCount = cumC[0];
        // update cumulative total
        for (int i = 0; i < capacity; i++) {
            cumC[i] = cumC[i] - removedCount;
        }
        if (capacity == 1) {
            cumC[0] = countAtNewSite;
        } else {
            final int[] cumCCopy = new int[capacity];
            System.arraycopy(cumC, 1, cumCCopy, 0, indexAtLast);
            cumC = cumCCopy;
            cumCCopy[indexAtLast] = cumCCopy[indexAtLast - 1] + countAtNewSite;
        }
    }

    private void fillUpCumulativeArray(final int countAtNewSite) {
        if (capacityMonitor == 0) {
            cumC[capacityMonitor] = countAtNewSite;
        } else {
            cumC[capacityMonitor] = cumC[capacityMonitor - 1] + countAtNewSite;
        }
    }

    /**
     * Returns true if the cumulativeSum array is full
     *
     * @return true if the cumulativeSum array is full
     */
    private boolean checkCumulativeCapacityReached() {
        return capacityMonitor == capacity;
    }

    public int getCumulativeCountTotal() {
        return cumC[capacity - 1];
    }

    public String toString() {
        final StringBuilder stringResult = new StringBuilder("[");
        for (int x = 0; x < capacityMonitor; x++) {
            stringResult.append('\t');
            stringResult.append(cumC[x]);
        }
        stringResult.append("\t]");
        return stringResult.toString();
    }


}
