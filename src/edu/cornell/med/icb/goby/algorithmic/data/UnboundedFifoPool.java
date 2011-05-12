/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.algorithmic.data;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;

/**
 * A first in first out pool backed by a fastutil ObjectArrayList.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2011
 *         Time: 1:24:38 PM
 */
public class UnboundedFifoPool<T> {
    ObjectArrayList<T> array;
    private int tailIndex = 0;
    private int headIndex = 0;
    int numElements = 0;

    /**
     * Create an UnboundedFifoPool with the specified initial capacity.
     *
     * @param capacity Number of elements the pool can initially hold. The capacity will increase as elements are added to the pool.
     */
    public UnboundedFifoPool(int capacity) {
        array = new ObjectArrayList<T>(capacity);
    }

    /**
     * Create an UnboundedFifoPool with default initial capacity.
     */
    public UnboundedFifoPool() {
        array = new ObjectArrayList<T>();
    }

    /**
     * Add an element to the pool.
     * @param element element to add.
     */
    public final void add(final T element) {
        //     System.out.printf("Adding %s pre: head-index=%d tail-index=%d %n", element, headIndex, tailIndex);

        ++numElements;
        array.ensureCapacity(numElements);
        array.size(Math.max(numElements, array.size()));
        array.set(tailIndex, element);
        advanceTailIndex();

        //  System.out.printf("Adding post: head-index=%d tail-index=%d array-size=%d %n", headIndex, tailIndex,
        //         array.size());
    }

    /**
     * Remove an element from the pool. The first element added will be removed.
     * @return
     */
    public final T remove() {
        //   System.out.printf("Removing pre: head-index=%d tail-index=%d %n", headIndex, tailIndex);
        if (isEmpty()) throw new IllegalStateException("Cannot remove element from empty pool");


        T element = array.get(headIndex);
        advanceHeadIndex();
        --numElements;

        if (isEmpty()) {
            headIndex = 0;
            tailIndex = 0;
        }
        //     System.out.printf("Removing post: head-index=%d tail-index=%d %n", headIndex, tailIndex);

        return element;

    }

    private void advanceHeadIndex() {

        ++headIndex;
        if (headIndex > numElements) {
            headIndex = 0;
        }
    }

    private void advanceTailIndex() {

        ++tailIndex;
        if (tailIndex > numElements) {
            tailIndex = 0;
        }
    }

    /**
     * Determine if the FIFO pool is empty.
     * @return True if the pool is empty, false otherwise.
     */
    public final boolean isEmpty() {
        return numElements == 0;
    }


}
