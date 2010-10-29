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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.counts.CountsReader;
import edu.cornell.med.icb.goby.counts.CountsWriter;
import edu.cornell.med.icb.util.RandomAdapter;
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
import java.util.Random;

/**
 * @author Fabien Campagne
 *         Date: Jun 5, 2009
 *         Time: 6:55:27 PM
 */
public class TestComputeStartCount {
    private static final Log LOG = LogFactory.getLog(TestAccumulate.class);
    private static final String BASE_TEST_DIR = "test-results/start-counts";

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

    @Test
    public void testComputeStarts() throws IOException {
        final ComputeStartCount computer = new ComputeStartCount(ComputeStartCount.POSITIVE_STRAND_ONLY);
        computer.populate(1, 10, true);
        computer.populate(1, 12, true);
        computer.populate(2, 3, true);
        computer.populate(2, 4, true);
        computer.populate(2, 12, true);

        final String filename = "test-results/start-counts/101.bin";
        final CountsWriter writer = new CountsWriter(new FileOutputStream(filename));
        computer.baseCount(writer);
        writer.close();

        final CountsReader reader = new CountsReader(new FileInputStream(filename));
        assertTrue(reader.hasNextPosition());
        assertEquals(0, reader.nextCountAtPosition());

        assertEquals(2, reader.nextCountAtPosition());
        assertEquals(3, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertFalse(reader.hasNextPosition());
        reader.close();
    }

    @Test
    public void testComputeStarts2() throws IOException {
        final ComputeStartCount computer = new ComputeStartCount(ComputeStartCount.POSITIVE_STRAND_ONLY);
        computer.populate(1, 10, true);
        computer.populate(1, 12, true);
        computer.populate(4, 30, true);
        computer.populate(4, 4, true);
        computer.populate(4, 12, true);
        computer.populate(10, 12, true);
        computer.populate(10, 12, true);
        computer.populate(11, 12, true);
        computer.populate(11, 12, true);
        computer.populate(12, 13, true);

        final String filename = "test-results/start-counts/101.bin";
        final CountsWriter writer = new CountsWriter(new FileOutputStream(filename));
        computer.baseCount(writer);
        writer.close();

        final CountsReader reader = new CountsReader(new FileInputStream(filename));
        assertTrue(reader.hasNextPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(0, reader.getPosition());
        assertEquals(2, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(3, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(9, reader.getPosition());
        assertEquals(2, reader.nextCountAtPosition());
        assertEquals(10, reader.getPosition());
        assertEquals(2, reader.nextCountAtPosition());
        assertEquals(11, reader.getPosition());
        assertEquals(1, reader.nextCountAtPosition());
        assertEquals(12, reader.getPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertFalse(reader.hasNextPosition());
        reader.close();
    }

    @Test
    public void testComputeStarts3() throws IOException {
        final ComputeStartCount computer = new ComputeStartCount(ComputeStartCount.POSITIVE_STRAND_ONLY);
        computer.populate(1, 10, true);
        computer.populate(1, 12, true);
        computer.populate(6, 3, true);
        computer.populate(6, 4, true);
        computer.populate(7, 8, true);
        computer.populate(7, 8, true);
        computer.populate(9, 12, true);
        computer.populate(9, 12, true);

        final String filename = "test-results/start-counts/101.bin";
        final CountsWriter writer = new CountsWriter(new FileOutputStream(filename));
        computer.baseCount(writer);
        writer.close();

        final CountsReader reader = new CountsReader(new FileInputStream(filename));
        assertTrue(reader.hasNextPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(2, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());

        assertEquals(2, reader.nextCountAtPosition());
        assertEquals(2, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertEquals(2, reader.nextCountAtPosition());
        assertEquals(0, reader.nextCountAtPosition());
        assertFalse(reader.hasNextPosition());
        reader.close();
    }

    @Test
    public void testComputeStarts4() throws IOException {
        final ComputeStartCount computer = new ComputeStartCount(ComputeStartCount.POSITIVE_STRAND_ONLY);
        for (int i = 0; i < 10000; i++) {
            if (i % 2 == 1) {
                computer.populate(i, 10001, true);
            }

        }


        final String filename = "test-results/start-counts/102.bin";
        final CountsWriter writer = new CountsWriter(new FileOutputStream(filename));
        computer.baseCount(writer);
        writer.close();

        final CountsReader reader = new CountsReader(new FileInputStream(filename));
        for (int i = 0; i < 10000; i++) {
            assertTrue(reader.hasNextPosition());

            final int count = reader.nextCountAtPosition();
            assertEquals((i % 2 == 1) ? 1 : 0, count);
            assertEquals(i, reader.getPosition());
        }

        assertTrue(reader.hasNextPosition());

        final int count = reader.nextCountAtPosition();
        assertEquals(0, count);
        assertFalse(reader.hasNextPosition());
        reader.close();
    }

    /**
     * @param lo lower limit of range
     * @param hi upper limit of range
     * @return a random integer in the range <STRONG>lo</STRONG>,
     *         <STRONG>lo</STRONG>+1, ... ,<STRONG>hi</STRONG>
     */
    private int chooseRandom(Random random, final int lo, final int hi) {
        final double r = random.nextDouble();
        int result = (int) ((long) lo + (long) ((1L + (long) hi - (long) lo) * r));
        assert result >= lo && result <= hi;
        return result;
    }

    @Test
    public void testComputeStarts5() throws IOException {
        initializeTestDirectory();
        final Random random = new Random();
        final ComputeStartCount computer = new ComputeStartCount(ComputeStartCount.POSITIVE_STRAND_ONLY);
        for (int i = 0; i < 100000; i++) {
            final int start = chooseRandom(random, 1, 10000);

            final int length = chooseRandom(random, 10, 100);
            computer.populate(start, start + length, true);

        }
        computer.populate(10010, 10020, true);
        computer.populate(10010, 10020, true);
        computer.populate(10011, 10020, true);

        final String filename = "test-results/start-counts/103.bin";
        final CountsWriter writer = new CountsWriter(new FileOutputStream(filename));
        computer.baseCount(writer);
        writer.close();

        final CountsReader reader = new CountsReader(new FileInputStream(filename));
        reader.skipTo(10003);


        assertEquals(2, reader.getCount());
        assertEquals(10010, reader.getPosition());
        reader.nextTransition();
        assertEquals(1, reader.getCount());
        assertEquals(10011, reader.getPosition());
        reader.close();
    }
}
