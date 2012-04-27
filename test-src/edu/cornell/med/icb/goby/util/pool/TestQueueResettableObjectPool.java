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

import it.unimi.dsi.fastutil.longs.LongArraySet;
import it.unimi.dsi.fastutil.longs.LongSet;
import org.junit.Test;

import java.util.Date;
import java.util.Random;

import static junit.framework.Assert.assertEquals;

/**
 * Test the queue object pool.
 */
public class TestQueueResettableObjectPool {
    /**
     * We use a seeded random number generator so we know in the first 100 there aren't repeats
     * which might cause the test to sometimes pass and sometimes fail.
     */
    private static final Random RANDOM = new Random(997);

    /**
     * Test the pool behaves as expected with regards to creation, pooling, and resetting().
     */
    @Test
    public void testPool() {
        final ResettableObjectPoolInterface<LocalResettable> pool = new QueueResettableObjectPool<LocalResettable>() {
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
        // Reset() is called upon return. It assumed if the object needs to be reset() upon creation,
        // the object will do it
        assertEquals("should be reset", 10, second.x);

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
        assertEquals("should be reset", 0, first.x);
        assertEquals("should be reset", 0, second.x);

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

    /**
     * Verify that the pool only creates new objects when expected.
     */
    @Test
    public void testAlwaysNewObjects() {
        final ResettableObjectPoolInterface<LocalResettable> pool = new QueueResettableObjectPool<LocalResettable>() {
            @Override
            public LocalResettable makeObject() {
                return new LocalResettable();
            }
        };

        final LongSet checkedOutRandoms = new LongArraySet();

        int sizeToExpect = 2;
        for (int i = 0; i < 50; i++) {
            final LocalResettable first = pool.borrowObject();  checkedOutRandoms.add(first.rand);
            final LocalResettable second = pool.borrowObject();  checkedOutRandoms.add(second.rand);
            if (RANDOM.nextBoolean()) {
                pool.returnObject(first);
                pool.returnObject(second);
            } else {
                if (RANDOM.nextBoolean()) {
                    pool.returnObject(second);
                    pool.returnObject(first);
                } else {
                    pool.returnObject(first);
                    pool.invalidateObject(second);
                    if (i != 49) {
                        sizeToExpect++;
                    }
                }
            }
        }
        assertEquals("Should be 100 objects created", sizeToExpect, checkedOutRandoms.size());
    }


    static class LocalResettable implements Resettable {
        int x = 10;
        long rand = RANDOM.nextLong();

        @Override
        public void reset() {
            x = 0;
        }
    }
}
