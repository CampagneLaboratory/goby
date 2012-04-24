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

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * Test the queue object pool.
 */
public class TestQueueObjectPool {
    @Test
    public void testPool() {
        final QueueObjectPool<LocalResettable> pool = new QueueObjectPool<LocalResettable>() {
            @Override
            public LocalResettable makeObject() {
                return new LocalResettable();
            }
        };

        LocalResettable first;
        assertEquals("Should be empty", 0, pool.getNumIdle());
        assertEquals("Should be empty", 0, pool.getNumActive());

        first = pool.borrowObject();
        first.x = 5;
        assertEquals("Should be empty", 0, pool.getNumIdle());
        assertEquals("Should be empty", 1, pool.getNumActive());

        pool.returnObject(first);
        assertEquals("Should be empty", 1, pool.getNumIdle());
        assertEquals("Should be empty", 0, pool.getNumActive());

        first = pool.borrowObject();
        LocalResettable second = pool.borrowObject();
        assertEquals("Should be empty", 0, pool.getNumIdle());
        assertEquals("Should be empty", 2, pool.getNumActive());
        assertEquals("should be reset", 0, first.x);
        assertEquals("should be reset", 0, second.x);

        pool.returnObject(second);
        assertEquals("Should be empty", 1, pool.getNumIdle());
        assertEquals("Should be empty", 1, pool.getNumActive());
        pool.returnObject(first);
        assertEquals("Should be empty", 2, pool.getNumIdle());
        assertEquals("Should be empty", 0, pool.getNumActive());


        first = pool.borrowObject();
        second = pool.borrowObject();
        assertEquals("Should be empty", 0, pool.getNumIdle());
        assertEquals("Should be empty", 2, pool.getNumActive());
        pool.invalidateObject(second);
        assertEquals("Should be empty", 0, pool.getNumIdle());
        assertEquals("Should be empty", 1, pool.getNumActive());
        pool.returnObject(first);
        assertEquals("Should be empty", 1, pool.getNumIdle());
        assertEquals("Should be empty", 0, pool.getNumActive());

        pool.clear() ;
        assertEquals("Should be empty", 0, pool.getNumIdle());
        assertEquals("Should be empty", 0, pool.getNumActive());

        first = pool.borrowObject();
        assertEquals("Should be empty", 0, pool.getNumIdle());
        assertEquals("Should be empty", 1, pool.getNumActive());
        pool.clear() ;
        assertEquals("Should be empty", 0, pool.getNumIdle());
        assertEquals("Should be empty", 1, pool.getNumActive());
        pool.returnObject(first);
        assertEquals("Should be empty", 1, pool.getNumIdle());
        assertEquals("Should be empty", 0, pool.getNumActive());

    }

    static class LocalResettable implements Resettable {
        int x = 10;
        @Override
        public void reset() {
            x = 0;
        }
    }
}
