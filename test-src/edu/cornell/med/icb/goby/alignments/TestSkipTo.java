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

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 6:15:52 PM
 */
public class TestSkipTo {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestSkipTo.class);
    private static final String BASE_TEST_DIR = "test-results/alignments-skip-to";
    private int numEntriesPerChunk=2;

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
    public void testFewSkips1() throws IOException {
        final String basename = "align-skip-to-1";
        final AlignmentWriter writer =
                new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);

        writer.setAlignmentEntry(0, 1, 12, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 2, 123, 30, false);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false);
        writer.appendEntry();
        writer.close();
        writer.printStats(System.out);

        final AlignmentReader reader =
                new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));

        final Alignments.AlignmentEntry a = reader.skipTo(0, 0);
        assertNotNull(a);
        assertEquals(1, a.getTargetIndex());
        assertEquals(12, a.getPosition());
        final Alignments.AlignmentEntry b = reader.skipTo(0, 0);
        assertEquals(1, b.getTargetIndex());
        assertEquals(13, b.getPosition());

        final Alignments.AlignmentEntry c = reader.skipTo(0, 0);
        assertEquals(1, c.getTargetIndex());
        assertEquals(13, c.getPosition());

        final Alignments.AlignmentEntry d = reader.skipTo(2, 300);
        assertEquals(2, d.getTargetIndex());
        assertEquals(300, d.getPosition());

        final Alignments.AlignmentEntry e = reader.skipTo(2, 300);
        assertEquals(2, e.getTargetIndex());
        assertEquals(300, e.getPosition());
        assertFalse(reader.hasNext());


    }

    @Test
    public void testFewSkips1_WithConcat() throws IOException {
        final String basename = "align-skip-to-1-concat";
        final AlignmentWriter writer =
                new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);

        writer.setAlignmentEntry(0, 1, 12, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 2, 123, 30, false);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false);
        writer.appendEntry();
        writer.close();
        writer.printStats(System.out);

        final ConcatSortedAlignmentReader reader = new ConcatSortedAlignmentReader(
                FilenameUtils.concat(BASE_TEST_DIR, basename),
                FilenameUtils.concat(BASE_TEST_DIR, basename));


        Alignments.AlignmentEntry a = reader.skipTo(0, 0);
        assertNotNull(a);
        assertEquals(1, a.getTargetIndex());
        assertEquals(12, a.getPosition());
        a = reader.skipTo(0, 0);
        assertNotNull(a);
        assertEquals(1, a.getTargetIndex());
        assertEquals(12, a.getPosition());

        final Alignments.AlignmentEntry b = reader.skipTo(0, 0);
        assertEquals(1, b.getTargetIndex());
        assertEquals(13, b.getPosition());

        Alignments.AlignmentEntry c = reader.next();
        assertEquals(1, c.getTargetIndex());
        assertEquals(13, c.getPosition());

        c = reader.next();
        assertEquals(1, c.getTargetIndex());
        assertEquals(13, c.getPosition());
        c = reader.next();
        assertEquals(1, c.getTargetIndex());
        assertEquals(13, c.getPosition());

        Alignments.AlignmentEntry d = reader.skipTo(2, 300);
        assertEquals(2, d.getTargetIndex());
        assertEquals(300, d.getPosition());

        d = reader.skipTo(2, 300);
        assertEquals(2, d.getTargetIndex());
        assertEquals(300, d.getPosition());

        Alignments.AlignmentEntry e = reader.skipTo(2, 300);
        assertNotNull(e);
        assertEquals(2, e.getTargetIndex());
        assertEquals(300, e.getPosition());
        e = reader.skipTo(2, 300);
        assertEquals(2, e.getTargetIndex());
        assertEquals(300, e.getPosition());
        assertFalse(reader.hasNext());


    }

    @Test
    public void testFewSkips2() throws IOException {
        final String basename = "align-skip-to-2";
        final AlignmentWriter writer =
                new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);

        writer.setAlignmentEntry(0, 1, 12, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 2, 123, 30, false);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false);
        writer.appendEntry();
        writer.close();
        writer.printStats(System.out);

        final AlignmentReader reader =
                new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));


        final Alignments.AlignmentEntry c = reader.skipTo(2, 0);
        assertEquals(2, c.getTargetIndex());
        assertEquals(123, c.getPosition());


        final Alignments.AlignmentEntry d = reader.skipTo(2, 300);
        assertEquals(2, d.getTargetIndex());
        assertEquals(300, d.getPosition());

        final Alignments.AlignmentEntry e = reader.skipTo(2, 300);
        assertEquals(2, e.getTargetIndex());
        assertEquals(300, e.getPosition());
        assertFalse(reader.hasNext());


    }

    @Test
    public void testEmptyAlignment() throws IOException {
        final String basename = "align-skip-to-3";
        final AlignmentWriter writer =
                new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);

        writer.close();
        writer.printStats(System.out);

        final AlignmentReader reader =
                new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));


        final Alignments.AlignmentEntry c = reader.skipTo(2, 0);
        assertNull(c);

        assertFalse(reader.hasNext());


    }


    @Test
    public void testFewSkips4() throws IOException {
        final String basename = "align-skip-to-2";
        final AlignmentWriter writer =
                new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, basename));
        writer.setNumAlignmentEntriesPerChunk(numEntriesPerChunk);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);

        writer.setAlignmentEntry(0, 1, 12, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 1, 13, 30, false);
        writer.appendEntry();

        writer.setAlignmentEntry(0, 2, 123, 30, false);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false);
        writer.appendEntry();
        writer.setAlignmentEntry(0, 2, 300, 30, false);
        writer.appendEntry();
        writer.close();
        writer.printStats(System.out);

        final AlignmentReader reader =
                new AlignmentReaderImpl(FilenameUtils.concat(BASE_TEST_DIR, basename));


        final Alignments.AlignmentEntry c = reader.skipTo(1, 13);
        assertEquals(1, c.getTargetIndex());
        assertEquals(13, c.getPosition());


        final Alignments.AlignmentEntry d = reader.next();
        assertEquals(1, d.getTargetIndex());
        assertEquals(13, d.getPosition());

        final Alignments.AlignmentEntry e = reader.next();
               assertEquals(2, e.getTargetIndex());
               assertEquals(123, e.getPosition());



    }

}
