/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.alignments;

import org.junit.Test;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.io.File;

import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntArrayList;

/**
 * @author Fabien Campagne
 *         Date: Jun 22, 2010
 *         Time: 11:26:20 AM
 */
public class TestConcatSortedAlignmentReader {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestReadWriteAlignments.class);
    private static final String BASE_TEST_DIR = "test-results/sort-concat";

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
    public void testSortConcat() throws IOException {
        ConcatSortedAlignmentReader concat = new ConcatSortedAlignmentReader(basename1, basename2, basename3);
        IntList sortedPositions = new IntArrayList();
        int[] expectedPositions = {1, 2, 3, 5, 6, 7, 8, 9, 10, 10, 12, 99};
        for (Alignments.AlignmentEntry entry : concat) {
            // System.out.println("entry.position(): "+entry.getPosition());
            sortedPositions.add(entry.getPosition());
        }
        assertEquals(expectedPositions.length, sortedPositions.size());
        for (int i = 0; i < expectedPositions.length; i++) {
            assertEquals(expectedPositions[i], sortedPositions.getInt(i));
        }
    }

    final static String basename1 = FilenameUtils.concat(BASE_TEST_DIR, "sort-concat-1");
    final static String basename2 = FilenameUtils.concat(BASE_TEST_DIR, "sort-concat-2");
    final static String basename3 = FilenameUtils.concat(BASE_TEST_DIR, "sort-concat-3");

    @Before
    public void setUp() throws IOException {

        final AlignmentWriter writer1 =
                new AlignmentWriter(basename1);
        writer1.setNumAlignmentEntriesPerChunk(1000);

        // we write this alignment sorted:
        writer1.setSorted(true);

        append(writer1, 1, 1);
        append(writer1, 1, 2);
        append(writer1, 1, 3);
        append(writer1, 1, 10);
        append(writer1, 1, 99);

        writer1.close();

        final AlignmentWriter writer2 =
                new AlignmentWriter(basename2);
        writer2.setNumAlignmentEntriesPerChunk(1000);

        // we write this alignment sorted:
        writer2.setSorted(true);

        append(writer2, 1, 5);
        append(writer2, 1, 8);
        append(writer2, 1, 9);
        append(writer2, 1, 12);

        writer2.close();

        final AlignmentWriter writer3 =
                new AlignmentWriter(basename3);
        writer3.setNumAlignmentEntriesPerChunk(1000);

        // we write this alignment sorted:
        writer3.setSorted(true);

        append(writer3, 1, 6);
        append(writer3, 1, 7);
        append(writer3, 1, 10);

        writer3.close();


    }

    private void append(AlignmentWriter writer, int referenceIndex, int position) throws IOException {
        writer.setAlignmentEntry(0, referenceIndex, position, 1, false);
        writer.appendEntry();
    }


}
