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
 * Facilitates adding new counts to the rightmost element in the array
 *
 * @Author: Nyasha Chambwe
 * Date: 3/1/12
 * Time: 4:14 PM
 */
public class SlidingCountArray {


    public void setCumC(int[] cumC) {
        this.cumC = cumC;
    }

    /**
     * Cumulative count array
     */
    int[] cumC;
    /*
    *Number of elements that are in the array
    * **/
    int capacityMonitor;

    /*
    * Point to the first logical element in the array
    * */
    int head;

    int tail;
    /*
    Maximum size of the cumulative count array
     */
    private final int capacity;

    public SlidingCountArray(final int n) {
        if (!(n > 0)) {
            throw new IllegalArgumentException("The declared size of the cumulative count array must be at least 1");
        }
        cumC = new int[n];
        capacity = n;
        head = 0;
        tail = 0;
    }

    /*
    * Update cumulative count array with a new observed count
    * whilst maintaining the size of the array
    * */
    public void addToRight(final int countAtNewSite) {
        if (!cumulativeCapacityReached()) {
            fillUpCumulativeArray(countAtNewSite);
            capacityMonitor++;
            tail = advanceToNextIndex(tail);
        } else {
        shift(countAtNewSite);
            if (head != tail) {
                head = advanceToNextIndex(head);
                tail = advanceToNextIndex(tail);
            } else {
               head = advanceToNextIndex(head);
            }
        }
    }

    /*
    * Advance pointers whilst staying in bounds of the array
    * */
    private int advanceToNextIndex(final int currentIndex) {
        int pointer = currentIndex + 1;
        if (pointer == capacity) {
            pointer = 0;
        }
        return pointer;
    }

    /*
    * Updates the cumulative count array by subtracting count of the leftmost element
    * and adding the new cumulative sum to the rightmost element
    * */
    private void shift(final int countAtNewSite) {
        final int removedCount = cumC[head];
        // subtract count being removed from each element
        for (int i = 0; i < capacity; i++) {
            cumC[i] = cumC[i] - removedCount;
        }
        // update cumulative total
        if (capacity == 1) {
            cumC[tail] = countAtNewSite;
        } else {
            if (head == 0 && tail == 0) {
                cumC[tail] = cumC[capacity - 1] + countAtNewSite;
            } else {
               cumC[head] = cumC[tail] + countAtNewSite;
            }
        }
    }

    private void fillUpCumulativeArray(final int countAtNewSite) {
        if (capacityMonitor == 0) {
            cumC[capacityMonitor] = countAtNewSite;
        } else {
            cumC[tail] = cumC[tail - 1] + countAtNewSite;
        }
    }

    /**
     * Checks if the cumulativeSum array is full
     *
     * @return true if the cumulativeSum array is full
     */
    private boolean cumulativeCapacityReached() {
        return capacityMonitor == capacity;
    }

    public String cumulativeArrayToString(final int[] result) {
        final StringBuilder outputResult = new StringBuilder("[");
        int pointer= head;
        for(int i=0; i< capacityMonitor; i++){
            outputResult.append("\t");
            outputResult.append(cumC[pointer]);
            pointer =  advanceToNextIndex(pointer);
      }
              outputResult.append("\t]");
        return outputResult.toString();
    }

    public String getCumulativeArrayAsString() {
        return cumulativeArrayToString(cumC);
    }

    /*
    * Returns the cumulative sum over an array of size n
    * */
    public int getCumulativeSum(){
        return cumC[tail];
    }
}
