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

package edu.cornell.med.icb.goby.util.pool;

import org.apache.commons.pool.BaseObjectPool;

import java.util.LinkedList;
import java.util.Queue;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Pool of Resettable objects. The queue/pool size is unbounded so this shouldn't be used
 * when lots and lots of items will simultaneously be borrowed. The idle queue size will never
 * be decreased unless clear() is called. t.reset() will be called when the item is borrowed.
 *
 * **  THIS CLASS IS NOT THREAD SAFE. **
 */
public abstract class QueueObjectPool<T extends Resettable> extends BaseObjectPool<T> {

    /**
     * The pool of objects waiting to be borrowed.
     */
    private final Queue<T> queue = new LinkedList<T>();

    /**
     * The number of borrowed, active objects.
     */
    private int numActive;

    /**
     * Abstract class for making objects for the pool.
     * @return the new object
     */
    public abstract T makeObject();

    /**
     * Borrow an object
     * @return borrowed object
     */
    @Override
    public T borrowObject() {
        final T t;
        if (queue.isEmpty()) {
            t = makeObject();
        } else {
            t = queue.poll();
        }
        numActive++;
        t.reset();
        return t;
    }

    /**
     * This does NOT check if the object has already been returned.
     * Do not return the same object twice!
     * @param t object to return
     */
    @Override
    public void returnObject(final T t) {
        numActive--;
        queue.add(t);
    }

    /**
     * Return an object but do not place it back into the pool. For some reason
     * the object has become bad.
     * @param t object being returned but not placed into the queue
     */
    @Override
    public void invalidateObject(final T t) {
        // Don't return the object to the queue. Just let it drop.
        numActive--;
    }

    /**
     * Number of idle objects (in the pool)
     * @return Number of idle objects
     */
    @Override
    public int getNumIdle() {
        return queue.size();
    }

    /**
     * Number of active objects (checked out of the pool, not yet returned).
     * @return Number of active objects
     */
    @Override
    public int getNumActive() {
        return numActive;
    }

    /**
     * Clear the queue.
     */
    @Override
    public void clear() {
        queue.clear();
    }
}
