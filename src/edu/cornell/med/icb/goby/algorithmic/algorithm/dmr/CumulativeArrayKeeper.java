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
 * Faciliates adding new counts to the rightmost element in the array without changing
 *
 *
 * @Author: Nyasha Chambwe
 * Date: 3/1/12
 * Time: 4:14 PM
 */
public class CumulativeArrayKeeper {


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
    private int capacity;

    public CumulativeArrayKeeper(final int n) {
        if(!(n > 0))
        throw new IllegalArgumentException("The declared size of the cumulative count array must be at least 1");
        cumC = new int[n];
        capacityMonitor = 0;
        indexAtLast = capacity - 1;
        capacity =n;
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
        int removedCount = cumC[0];
        // update cumulative total
        for(int i=0; i < capacity; i++){
            cumC[i]=cumC[i]-removedCount;
        }
        if(capacity==1){
            cumC[0]= countAtNewSite;
        }else{
        int[] cumCCopy= new int[capacity];
        System.arraycopy(cumC, 1, cumCCopy, 0, indexAtLast);
        cumC= cumCCopy;
        cumCCopy[indexAtLast]=cumCCopy[indexAtLast-1]+countAtNewSite;
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
     * @return 
     */
    private boolean checkCumulativeCapacityReached() {
        if (capacityMonitor == capacity) return true;
        else return false;
    }

    public int getCumulativeCountTotal() {
        return cumC[capacity - 1];
    }

    public String cumulativeArrayToString(int[] result) {
        StringBuilder stringresult = new StringBuilder("[");
        int max = getCurrentCapacity();
        for (int x = 0; x < max; x++) {
            stringresult.append("\t");
            stringresult.append(cumC[x]);
        }
        stringresult.append("\t]");
        return stringresult.toString();
    }

    private int getCurrentCapacity() {
        return capacityMonitor;
    }

    public String getCumulativeArray() {
        return cumulativeArrayToString(this.cumC);
    }
}
