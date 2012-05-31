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
 * Largely based off of commons ObjectPool, but simplified, without exceptions, etc. Also, the objects
 * stored in pools based from this interface must be Resettable.
 */

public interface ResettableObjectPoolInterface<T extends Resettable> {

    /**
     * Abstract class for making objects for the pool / to be borrowed.
     * @return the new object
     */
    public T makeObject();

    /**
     * Borrow an object from the pool.
     * @return borrowed object
     */
    public T borrowObject();

    /**
     * Return an object to the pool that was borrowed.
     * @param t object to return
     */
    public void returnObject(final T t);

    /**
     * Return an object but do not place it back into the pool. For some reason
     * the object has become bad.
     * @param t object being returned but not placed into the queue
     */
    public void invalidateObject(final T t);

    /**
     * Number of idle objects (in the pool).
     * @return Number of idle objects
     */
    public int getNumIdle();

    /**
     * Number of active objects (checked out of the pool, not yet returned).
     * @return Number of active objects
     */
    public int getNumActive();

    /**
     * Clear objects stored in the pool waited to be borrowed.
     */
    public void clear();
}
