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

package edu.cornell.med.icb.goby.alignments.perms;

import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.apache.commons.io.FileUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 3/9/12
 *         Time: 2:50 PM
 */
public class TestPermutationWriter {
    private static final Log LOG = LogFactory.getLog(TestPermutationWriter.class);
    private static final String BASE_TEST_DIR = "test-results/permutations";

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
        // FileUtils.forceDeleteOnExit(new File(BASE_TEST_DIR));
    }

    @Test
    public void testBreakpoint() throws Exception {
        final PermutationWriter writer = new PermutationWriter("test-results/permutations/breakpoint-1");
        IntArrayList list = IntArrayList.wrap(new int[]{0, 1, 2, 4, 5});

        assertEquals(3, writer.getBreakPoint(0, list));
    }

    @Test
    public void testBreakpoint2() throws Exception {
        final PermutationWriter writer = new PermutationWriter("test-results/permutations/breakpoint-2");
        IntArrayList list = IntArrayList.wrap(new int[]{0, 1, 2, 3, 4, 5});

        assertEquals(3, writer.getBreakPoint(1, list, 2));
    }

    @Test
    public void testBreakpoint3() throws Exception {
        final PermutationWriter writer = new PermutationWriter("test-results/permutations/breakpoint-3");
        IntArrayList list = IntArrayList.wrap(new int[]{0, 1, 2, 3, 4, 5});

        assertEquals(6, writer.getBreakPoint(4, list, 20));
    }

    @Test
    public void testBreakpoint4() throws Exception {
        final PermutationWriter writer = new PermutationWriter("test-results/permutations/breakpoint-3");
        IntArrayList list = IntArrayList.wrap(new int[]{0, 1, 2, 3, 4, 5});

        assertEquals(6, writer.getBreakPoint(6, list, 20));
    }

    @Test
    public void testAppend() throws Exception {
        final String basename = "test-results/permutations/append-1";
        final PermutationWriter writer = new PermutationWriter(basename);
        Int2IntMap map = new Int2IntArrayMap();
        for (int queryIndex = 0; queryIndex < 100; queryIndex++) {
            final int smallIndex = queryIndex + 1;
            map.put(queryIndex, smallIndex);
        }
        for (int i = 300; i < 310; i++) {
            map.put(i, i + 1);
        }
        writer.append(map);
        writer.close();

        final PermutationReaderInterface reader = new PermutationReader(basename);
        for (int i = 1; i < 101; i++) {
            int expectedQueryIndex = i - 1;
            final int queryIndex = reader.getQueryIndex(i);
            assertEquals("got wrong answer for smallIndex=" + i, expectedQueryIndex, queryIndex);
        }
        for (int i = 301; i < 311; i++) {
            int expectedQueryIndex = i - 1;
            final int queryIndex = reader.getQueryIndex(i);
            assertEquals("got wrong answer for smallIndex=" + i, expectedQueryIndex, queryIndex);
        }

        assertEquals("got wrong answer", -1, reader.getQueryIndex(200));
        assertEquals("got wrong answer", -1, reader.getQueryIndex(250));
        assertEquals("got wrong answer", -1, reader.getQueryIndex(0));

    }
}
