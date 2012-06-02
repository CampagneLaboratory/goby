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

import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.longs.LongArraySet;
import it.unimi.dsi.fastutil.longs.LongList;
import it.unimi.dsi.fastutil.longs.LongSet;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Random;

import static junit.framework.Assert.assertEquals;

/**
 * Test the queue object pool.
 */
public class TestNullResettableObjectPool {
    /**
     * We use a seeded random number generator so we know in the first 100 there aren't repeats
     * which might cause the test to sometimes pass and sometimes fail.
     */
    private static final Random RANDOM = new Random(997);

    /**
     * Verify that we create always new objects when we check out and return
     */
    @Test
    public void testAlwaysNewObjects() {
        final ResettableObjectPoolInterface<LocalResettable> pool = new NullResettableObjectPool<LocalResettable>() {
            @Override
            public LocalResettable makeObject() {
                return new LocalResettable();
            }
        };

        final LongSet checkedOutRandoms = new LongArraySet();

        LocalResettable first, second;
        for (int i = 0; i < 50; i++) {
            first = pool.borrowObject();  checkedOutRandoms.add(first.rand);
            second = pool.borrowObject();  checkedOutRandoms.add(second.rand);
            assertEquals("Should be empty", 0, pool.getNumIdle());
            assertEquals("Should be empty", 0, pool.getNumActive());
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
                }
            }
            assertEquals("Should be empty", 0, pool.getNumIdle());
            assertEquals("Should be empty", 0, pool.getNumActive());
        }
        assertEquals("Should be 100 objects created", 100, checkedOutRandoms.size());
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
