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
 */
public abstract class QueueObjectPool<T extends Resettable> extends BaseObjectPool<T> {

    private final Queue<T> queue = new LinkedList<T>();

    private final AtomicInteger numActive = new AtomicInteger(0);

    public abstract T makeObject();

    @Override
    public synchronized T borrowObject() {
        final T t;
        if (queue.isEmpty()) {
            t = makeObject();
        } else {
            t = queue.poll();
        }
        numActive.incrementAndGet();
        t.reset();
        return t;
    }

    @Override
    public synchronized void returnObject(final T t) {
        if (!queue.contains(t)) {
            numActive.decrementAndGet();
            queue.add(t);
        }
    }

    @Override
    public synchronized void invalidateObject(final T t) {
        // Don't return the object to the queue. Just let it drop.
        numActive.decrementAndGet();
    }

    @Override
    public int getNumIdle() {
        return queue.size();
    }

    @Override
    public int getNumActive() {
        return numActive.get();
    }

    @Override
    public synchronized void clear() {
        queue.clear();
    }
}
