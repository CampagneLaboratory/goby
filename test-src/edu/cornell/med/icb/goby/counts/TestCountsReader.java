/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.counts;

import edu.cornell.med.icb.goby.algorithmic.algorithm.TestAccumulate;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.io.FastByteArrayInputStream;
import org.apache.commons.io.FileUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;
import java.util.Random;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author campagne
 *         Date: 6/15/11
 *         Time: 6:53 PM
 */
public class TestCountsReader {
    private static final Log LOG = LogFactory.getLog(TestAccumulate.class);
    private static final String BASE_TEST_DIR = "test-results/counts/";
    private CountsReader reader;

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

            FileUtils.forceDeleteOnExit(new File(BASE_TEST_DIR));
        }
    }

    @Before
    public void createCounts() throws IOException {
        ByteArrayOutputStream stream = new ByteArrayOutputStream(100000);
        ByteArrayOutputStream indexByteArrayOutputStream = new ByteArrayOutputStream(100000);
        DataOutput indexOutput = new DataOutputStream(indexByteArrayOutputStream);
        CountsWriter writer = new CountsWriter(stream);
        writer.appendCount(0, 10); // before- Position=0  count-after=0
        writer.appendCount(1, 10); // before- Position=10
        writer.appendCount(2, 10); // before- Position=20
        writer.appendCount(3, 10); // before- Position=30    count-after=3
        writer.appendCount(10, 10); // before- Position=40   count-after=10
        writer.appendCount(11, 10); // before- Position=50   count-after=11
        writer.appendCount(9, 10);  // before- Position=60
        writer.appendCount(7, 10);  // before- Position=70
        writer.close();

        byte[] bytes = stream.toByteArray();
        CountIndexBuilder indexBuilder = new CountIndexBuilder(1);
        indexBuilder.buildIndex(bytes, indexOutput);
        byte[] bytes1 = indexByteArrayOutputStream.toByteArray();

        reader = new CountsReader(new FastByteArrayInputStream(bytes),
                new DataInputStream(new ByteArrayInputStream(bytes1)));


    }


    @Test
    public void testSkipToBetweenIndexPositions() throws IOException {
        reader.skipTo(35);

        assertEquals(40, reader.getPosition());
        assertEquals(10, reader.getCount());

        // skip to a past position:
        reader.skipTo(35);

        assertEquals(50, reader.getPosition());
        assertEquals(11, reader.getCount());
        // skip to exactly position 60, in the index and in the future:


    }


    @Test
    public void testSkipToAtIndexPositions() throws IOException {

        reader.skipTo(60);

        assertEquals(60, reader.getPosition());
        assertEquals(9, reader.getCount());
        reader.nextTransition();
        assertEquals(70, reader.getPosition());
        assertEquals(7, reader.getCount());
    }

    @Test
    public void testSkipToBoth() throws IOException {

        reader.skipTo(35);

        assertEquals(40, reader.getPosition());
        assertEquals(10, reader.getCount());

        reader.skipTo(60);

        assertEquals(60, reader.getPosition());
        assertEquals(9, reader.getCount());
    }


    @Test
    public void testReposition() throws IOException {


        reader.reposition(60);

        assertEquals(60, reader.getPosition());
        assertEquals(9, reader.getCount());

        reader.reposition(30);

        assertEquals(30, reader.getPosition());
        assertEquals(3, reader.getCount());
        reader.reposition(70);

        assertEquals(70, reader.getPosition());
        assertEquals(7, reader.getCount());

        reader.reposition(15);

        assertEquals(20, reader.getPosition());
        assertEquals(2, reader.getCount());
    }

    @Test
    public void testNegativeSkipTo() throws IOException {


        final CountsArchiveReader archiveReader = new CountsArchiveReader("test-data/counts/GDFQPGI-pickrellNA18486_yale");
        final CountsReader reader1 = archiveReader.getCountReader("1");
        reader1.reposition(77120584);
        reader1.skipTo(77120584);

        // count should be 0
    }

    @Test
    public void testNotNegativeAfterSkipTo() throws IOException {

        final CountsArchiveReader archiveReader = new CountsArchiveReader("test-data/counts/GDFQPGI-pickrellNA18486_yale");
        final CountsReader reader1 = archiveReader.getCountReader("1");
        int maxPos = -1;
        Int2IntMap positionToCount = new Int2IntOpenHashMap();
        Int2IntMap positionToLength = new Int2IntOpenHashMap();
        int maxTransitions=10000;
        int counter=0;
        int reader1Position = reader1.getPosition();
        while (reader1.hasNextTransition()) {
            reader1.nextTransition();
            maxPos = Math.max(maxPos, reader1Position);
            int position = reader1.getPosition();
            int count = reader1.getCount();
            int length = reader1.getLength();
            for (int j = 0; j < length; j++) {
                positionToCount.put(position+ j, count );
            }

            positionToLength.put(position, length);
            if (counter++>maxTransitions) break;
        }
        System.out.printf("maxPos=%,d%n", maxPos);
        Random random = new Random();
        int lastSkipToPos = -1;
        try {
            for (int i = 1; i < maxPos; i++) {


                    int skipToPos = (int) (random.nextDouble() * (maxPos - 1));

                    lastSkipToPos = skipToPos;

                    reader1.reposition(skipToPos);
                    if (positionToCount.containsKey(reader1Position) && reader1.getCount() != positionToCount.get(reader1Position)) {
                        System.out.printf("error, count differ at position %,d. Count was %d, should have been %d %n",
                                reader1Position,
                                reader1.getCount(),
                                positionToCount.get(reader1Position));
                    }
                    reader1.skipTo(skipToPos);

                    if (reader1.hasNextTransition()) {
                        reader1.nextTransition();
                        lastSkipToPos = skipToPos;
                        assertTrue("count must be positive. ", reader1.getCount() >= 0);
                        assertTrue("reader position must be greater or equal to skipToPos. ", reader1Position >= skipToPos);
                    }
                }

        } finally {
            System.out.printf("last skipToPosition=%,d %n", lastSkipToPos);
        }
    }

}
