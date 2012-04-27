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

/**
 * This class implements the ResettableObjectPoolInterface, but when you borrow, you always create a new object and
 * when you return the object, reset is called on the object (to free resources, etc.) but the object
 * won't be pooled.
 */
public abstract class NullResettableObjectPool<T extends Resettable> implements ResettableObjectPoolInterface<T> {

    /**
     * Abstract class for making objects for the pool.
     * @return the new object
     */
    @Override
    public abstract T makeObject();

    /**
     * Borrow an object from the pool.
     * Since this implementation doesn't keep a pool, this always returns a new object.
     * @return borrowed object
     */
    @Override
    public T borrowObject() {
        return makeObject();
    }

    /**
     * Return an object to the pool that was borrowed.
     * reset() will immediately be called on the object.
     * Since this implementation doesn't keep a pool, this doesn't actually return the object to the pool.
     * @param t object to return
     */
    @Override
    public void returnObject(final T t) {
        t.reset();
    }

    /**
     * Return an object but do not place it back into the pool as for some reason the object has become bad.
     * reset() will immediately be called on the object.
     * Since this implementation doesn't keep a pool, this doesn't actually return the object to the pool.
     * @param t object being returned but not placed into the queue
     */
    @Override
    public void invalidateObject(final T t) {
        t.reset();
    }

    /**
     * Number of idle objects (in the pool).
     * Since this implementation doesn't keep a pool, always returns 0.
     * @return Number of idle objects
     */
    @Override
    public int getNumIdle() {
        return 0;
    }

    /**
     * Number of active objects (checked out of the pool, not yet returned).
     * Since this implementation doesn't keep a pool, always returns 0.
     * @return Number of active objects
     */
    @Override
    public int getNumActive() {
        return 0;
    }

    /**
     * Clear objects stored in the pool waited to be borrowed.
     * Since this implementation doesn't keep a pool, this does nothing.
     */
    @Override
    public void clear() {
    }
}
