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

import java.util.LinkedList;
import java.util.Queue;

/**
 * Pool of Resettable objects. The queue/pool size is unbounded so this shouldn't be used
 * when lots and lots of items will simultaneously be borrowed. The idle queue size will never
 * be decreased unless clear() is called. t.reset() at the point the object is returned/invalidated.
 *
 * **  THIS CLASS IS NOT THREAD SAFE. **
 */
public abstract class QueueResettableObjectPool<T extends Resettable> implements ResettableObjectPoolInterface<T> {

    /**
     * The pool of objects waiting to be borrowed.
     */
    private final Queue<T> queue = new LinkedList<T>();

    /**
     * Abstract class for making objects for the pool / to be borrowed.
     */
    private int numActive;

    /**
     * Abstract class for making objects for the pool.
     * @return the new object
     */
    @Override
    public abstract T makeObject();

    /**
     * Borrow an object from the pool.
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
        return t;
    }

    /**
     * Return an object to the pool that was borrowed.
     * reset() will immediately be called on the object.
     * This does NOT check if the object has already been returned, so  do not return the same object twice.
     * @param t object to return
     */
    @Override
    public void returnObject(final T t) {
        t.reset();
        queue.add(t);
        numActive--;
    }

    /**
     * Return an object but do not place it back into the pool as for some reason the object has become bad.
     * reset() will immediately be called on the object.
     * @param t object being returned but not placed into the queue
     */
    @Override
    public void invalidateObject(final T t) {
        // Don't return the object to the queue. Just let it drop.
        t.reset();
        numActive--;
    }

    /**
     * Number of idle objects (in the pool).
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
     * Clear objects stored in the pool waited to be borrowed.
     */
    @Override
    public void clear() {
        queue.clear();
    }
}
