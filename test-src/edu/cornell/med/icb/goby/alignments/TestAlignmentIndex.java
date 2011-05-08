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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.objects.ObjectList;
import static junit.framework.Assert.assertEquals;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: May 6, 2011
 *         Time: 4:48:24 PM
 */
public class TestAlignmentIndex {
    /**
     * This test triggers an issue with Goby 1.9.5 indices. Using pre 1.9.6 indices, the skipTo method could fail
     * to return alignments if they occured at a position after the length of the next sequence. This test will fail
     * with any version of Goby prior to 1.9.6, and will success with Goby 1.9.6. Old alignments can be upgraded to
     * the index structures used by Goby 1.9.6 with the goby upgrade mode.
     *
     * @throws IOException
     */
    @Test
    public void test_1_9_5_IndexIssue() throws IOException {

        int[] targetLengths = new int[]{100, 50, 20, 10, 5};
        final String basename1 = FilenameUtils.concat(BASE_TEST_DIR, "align-index-error-1");
        final AlignmentWriter writer =
                new AlignmentWriter(basename1);
        writer.setTargetLengths(targetLengths);
        writer.setNumAlignmentEntriesPerChunk(1);
        writer.setSorted(true);
        int queryIndex = 0;

        writer.setAlignmentEntry(queryIndex++, 0, 99, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(queryIndex++, 1, 0, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(queryIndex++, 1, 1, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.close();

        AlignmentReaderImpl reader = new AlignmentReaderImpl(basename1);
        reader.reposition(0, 99);
        Alignments.AlignmentEntry entry = reader.next();
        assertEquals(0, entry.getTargetIndex());
        assertEquals(99, entry.getPosition());
        entry = reader.next();
        assertEquals(1, entry.getTargetIndex());
        assertEquals(0, entry.getPosition());
        entry = reader.next();
        assertEquals(1, entry.getTargetIndex());
        assertEquals(1, entry.getPosition());
        assertFalse(reader.hasNext());
        // Now check that the locations were stored in the index and can be decoded correctly:
        ObjectList<ReferenceLocation> locations = reader.getLocations(1);
        assertEquals(0, locations.get(0).targetIndex);
        assertEquals(99, locations.get(0).position);
        assertEquals(1, locations.get(1).targetIndex);
        assertEquals(0, locations.get(1).position);
        assertEquals(1, locations.get(2).targetIndex);
        assertEquals(1, locations.get(2).position);

        assertEquals("with modulo=1, must recover three locations.", 3, locations.size());
    }

    @Test
    public void test_1_9_5_IndexIssueExtraTests() throws IOException {

        int[] targetLengths = new int[]{100, 50, 20, 10, 5};
        final String basename1 = FilenameUtils.concat(BASE_TEST_DIR, "align-index-error-2");
        final AlignmentWriter writer =
                new AlignmentWriter(basename1);
        writer.setTargetLengths(targetLengths);
        writer.setNumAlignmentEntriesPerChunk(1);
        writer.setSorted(true);
        int queryIndex = 0;

        writer.setAlignmentEntry(queryIndex++, 0, 99, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(queryIndex++, 1, 0, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(queryIndex++, 1, 1, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.close();

        AlignmentReaderImpl reader = new AlignmentReaderImpl(basename1);
        reader.reposition(1, 0);    // will be (0,99) because reposition goes one chunk before. The previous chunk
        // will have the indexed entry as first entry, plus some other entry, possibly including an entry with
        // position equal to the reposition argument. Repositioning to the chunk before garantees that skipTo will
        // not miss these entries, and appropriately filter the beginning of the chunk with position before (1,0).

        Alignments.AlignmentEntry entry = reader.next();
        assertEquals(0, entry.getTargetIndex());
        assertEquals(99, entry.getPosition());
        entry = reader.next();
        assertEquals(1, entry.getTargetIndex());
        assertEquals(0, entry.getPosition());
        entry = reader.next();
        assertEquals(1, entry.getTargetIndex());
        assertEquals(1, entry.getPosition());
        assertFalse(reader.hasNext());
        // Now check that the locations were stored in the index and can be decoded correctly:
        ObjectList<ReferenceLocation> locations = reader.getLocations(1);
        assertEquals(0, locations.get(0).targetIndex);
        assertEquals(99, locations.get(0).position);
        assertEquals(1, locations.get(1).targetIndex);
        assertEquals(0, locations.get(1).position);
        assertEquals(1, locations.get(2).targetIndex);
        assertEquals(1, locations.get(2).position);

        assertEquals("with modulo=1, must recover three locations.", 3, locations.size());

        reader = new AlignmentReaderImpl(basename1);

        entry =  reader.skipTo(1, 1); 
        assertEquals(1, entry.getTargetIndex());
        assertEquals(1, entry.getPosition());
        assertFalse(reader.hasNext());
    }

    @Test
    public void test_1_9_5_IndexIssueWithSkipTo() throws IOException {

        int[] targetLengths = new int[]{100, 50, 20, 10, 5};
        final String basename1 = FilenameUtils.concat(BASE_TEST_DIR, "align-index-error-2");
        final AlignmentWriter writer =
                new AlignmentWriter(basename1);
        writer.setTargetLengths(targetLengths);
        writer.setNumAlignmentEntriesPerChunk(1);
        writer.setSorted(true);
        int queryIndex = 0;

        writer.setAlignmentEntry(queryIndex++, 0, 99, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(queryIndex++, 1, 0, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(queryIndex++, 1, 1, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.close();

        AlignmentReaderImpl reader = new AlignmentReaderImpl(basename1);
        reader.readHeader();
        Alignments.AlignmentEntry entry = reader.skipTo(1, 0);

        assertEquals(1, entry.getTargetIndex());
        assertEquals(0, entry.getPosition());
        entry = reader.next();
        assertEquals(1, entry.getTargetIndex());
        assertEquals(1, entry.getPosition());
        assertFalse(reader.hasNext());
        // Now check that the locations were stored in the index and can be decoded correctly:
        ObjectList<ReferenceLocation> locations = reader.getLocations(1);
        assertEquals(0, locations.get(0).targetIndex);
        assertEquals(99, locations.get(0).position);
        assertEquals(1, locations.get(1).targetIndex);
        assertEquals(0, locations.get(1).position);
        assertEquals(1, locations.get(2).targetIndex);
        assertEquals(1, locations.get(2).position);

        assertEquals("with modulo=1, must recover three locations.", 3, locations.size());
    }
    @Test
    public void testReposition() throws IOException {
         int[] targetLengths = new int[]{100, 50, 20, 10, 5};
        final String basename1 = FilenameUtils.concat(BASE_TEST_DIR, "align-index-error-2");
        final AlignmentWriter writer =
                new AlignmentWriter(basename1);
        writer.setTargetLengths(targetLengths);
        writer.setNumAlignmentEntriesPerChunk(1);
        writer.setSorted(true);
        int queryIndex = 0;

        writer.setAlignmentEntry(queryIndex++, 0, 99, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(queryIndex++, 1, 0, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.setAlignmentEntry(queryIndex++, 1, 1, 30, false, constantQueryLength);
        writer.appendEntry();
        writer.close();

        AlignmentReaderImpl reader = new AlignmentReaderImpl(basename1);
        reader.readHeader();
        assertNotNull(reader.skipTo(1,0));
        // now reposition to an earlier position (since 1.9.6):
        reader.reposition(0,99);
        assertNotNull(reader.skipTo(1,0));
        assertNotNull(reader.skipTo(1,1));
             reader.reposition(0,99);
         assertNotNull(reader.skipTo(1,0));
        assertNotNull(reader.skipTo(1,1));

    }
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestReadWriteAlignments.class);
    private static final String BASE_TEST_DIR = "test-results/alignment-index";
    private int constantQueryLength = 40;

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

}
