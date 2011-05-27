/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

package edu.cornell.med.icb.goby.counts;

import edu.cornell.med.icb.goby.algorithmic.algorithm.ComputeStartCount;
import org.apache.commons.io.FileUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Jun 13, 2009
 *         Time: 2:50:19 PM
 */
public class TestAnyTransitionCountsIterator {
    private static final Log LOG = LogFactory.getLog(TestAnyTransitionCountsIterator.class);
    private static final String BASE_TEST_DIR = "test-results/any-transition";

    @BeforeClass
    public static void initializeTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Creating base test directory: " + BASE_TEST_DIR);
        }
        FileUtils.forceMkdir(new File(BASE_TEST_DIR));
    }

    @AfterClass
    public static void cleanupTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Deleting base test directory: " + BASE_TEST_DIR);
        }
        FileUtils.forceDeleteOnExit(new File(BASE_TEST_DIR));
    }
     /*
    @Test
    public void testSimple1() throws IOException {
        final int[] starts1 = {
                1,
                2, 2, 2,
                4, 4,
                // no transition at position 7 for this reader.
                9,
        };

        final int[] starts2 = {
                1,
                2, 2, 2, 2,
                // no transition at position 4 for this reader..
                7,
                9,
                10
        };

        final String filename1 = "test-results/any-transition/" + "counts-1.bin";
        makeStartCounts(starts1, filename1);
        final String filename2 = "test-results/any-transition/" + "counts-2.bin";
        makeStartCounts(starts2, filename2);

        AnyTransitionCountsIterator iterator = null;
        try {
            final CountsReader reader1 = new CountsReader(new FileInputStream(filename1));
            final CountsReader reader2 = new CountsReader(new FileInputStream(filename2));
            iterator = new AnyTransitionCountsIterator(reader1, reader2);

            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("first position must be zero", 0, iterator.getPosition());
            assertEquals("reader 0 must have zero count", 0, iterator.getCount(0));
            assertEquals("reader 1 must have zero count", 0, iterator.getCount(1));

            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("second position must be 1", 1, iterator.getPosition());
            assertEquals("reader 0 must have count=1", 1, iterator.getCount(0));
            assertEquals("reader 1 must have count=1", 1, iterator.getCount(1));

            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("next position must be 2", 2, iterator.getPosition());
            assertEquals("reader 0 must have count=3", 3, iterator.getCount(0));
            assertEquals("reader 1 must have count=4", 4, iterator.getCount(1));

            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("next position must be 3", 3, iterator.getPosition());
            assertEquals("reader 0 must have count=0", 0, iterator.getCount(0));
            assertEquals("reader 1 must have count=0", 0, iterator.getCount(1));

            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("next position must be 4", 4, iterator.getPosition());
            assertEquals("reader 0 must have count=2", 2, iterator.getCount(0));
            assertEquals("reader 1 must have count=0", 0, iterator.getCount(1));  // reader 1 still has count=0 at position 4
            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("next position must be 5", 5, iterator.getPosition());
            assertEquals("reader 0 must have count=0", 0, iterator.getCount(0));  // reader 0 is back to zero
            assertEquals("reader 1 must have count=1", 0, iterator.getCount(1));

            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("next position must be 7", 7, iterator.getPosition());
            assertEquals("reader 0 must have count=0", 0, iterator.getCount(0));  // reader 0 still has count=0 at position 7
            assertEquals("reader 1 must have count=1", 1, iterator.getCount(1));

            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("next position must be 8", 8, iterator.getPosition());
            assertEquals("reader 0 must have count=0", 0, iterator.getCount(0));
            assertEquals("reader 1 must have count=1", 0, iterator.getCount(1));  // reader 0 is back to zero

            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("next position must be 9", 9, iterator.getPosition());
            assertEquals("reader 0 must have count=0", 1, iterator.getCount(0));
            assertEquals("reader 1 must have count=1", 1, iterator.getCount(1));

            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("next position must be 10", 10, iterator.getPosition());
            assertEquals("reader 0 must have count=0", 0, iterator.getCount(0));
            assertEquals("reader 1 must have count=1", 1, iterator.getCount(1));

            assertTrue(iterator.hasNextTransition());
            iterator.nextTransition();
            assertEquals("next position must be 11", 11, iterator.getPosition());
            assertEquals("reader 0 must have count=0", 0, iterator.getCount(0));
            assertEquals("reader 1 must have count=1", 0, iterator.getCount(1));

            assertFalse(iterator.hasNextTransition());
        } finally {
            if (iterator != null) {
                iterator.close();
            }
        }

    }
      */
    private void makeStartCounts(final int[] starts, final String filename) throws IOException {
        final ComputeStartCount computer = new ComputeStartCount(ComputeStartCount.POSITIVE_STRAND_ONLY);
        for (final int start : starts) {
            computer.populate(start, 10, true);
        }

        final CountsWriter writer = new CountsWriter(new FileOutputStream(filename));
        computer.baseCount(writer);
        writer.close();
    }
}
