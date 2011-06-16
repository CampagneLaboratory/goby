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

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.*;

/**
 * @author Fabien Campagne
 *         Date: May 14, 2009
 *         Time: 11:58:32 AM
 */
public class TestCountsArchive {
    private static final Log LOG = LogFactory.getLog(TestCountsArchive.class);
    private static final String BASE_TEST_DIR = "test-results/multicount";

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
    public void testWriter() throws IOException {
        final String basename = FilenameUtils.concat(BASE_TEST_DIR, "101.bin");
        final CountsArchiveWriter writer = new CountsArchiveWriter(basename);
        CountsWriter cw = writer.newCountWriter(0, "count-0");
        cw.appendCount(0, 10000);
        cw.appendCount(10, 10000);
        cw.appendCount(20, 10000);
        cw.close();
        writer.returnWriter(cw);

        // read counts 1
        cw = writer.newCountWriter(1, "count-1");
        cw.appendCount(0, 20000);
        cw.appendCount(10, 20000);
        cw.appendCount(20, 20000);
        cw.appendCount(30, 20000);
        cw.close();
        writer.returnWriter(cw);
        writer.close();

        final CountsArchiveReader reader = new CountsArchiveReader(basename);
        assertTrue("Since Goby 1.9.7 reader must have stats.", reader.isStatsParsed());
        CountsReader cr = reader.getCountReader("count-1");
        assertTrue(cr.hasNextTransition());
        cr.nextTransition();
        assertTrue(cr.hasNextTransition());
        cr.nextTransition();
        assertTrue(cr.hasNextTransition());
        cr.nextTransition();
        assertEquals(20, cr.getCount());
        assertTrue(cr.hasNextTransition());

        cr.nextTransition();
        assertEquals(30, cr.getCount());
        assertFalse(cr.hasNextTransition());
        cr.close();
        // read counts 0

        cr = reader.getCountReader("count-0");
        assertTrue(cr.hasNextTransition());
        cr.nextTransition();
        assertTrue(cr.hasNextTransition());
        cr.nextTransition();
        assertTrue(cr.hasNextTransition());
        cr.nextTransition();
        assertEquals(20, cr.getCount());

        assertFalse(cr.hasNextTransition());
        cr.close();

        // and again counts-1
        cr = reader.getCountReader("count-1");
        assertTrue(cr.hasNextTransition());
        cr.nextTransition();
        assertTrue(cr.hasNextTransition());
        cr.nextTransition();
        assertTrue(cr.hasNextTransition());
        cr.nextTransition();
        assertEquals(20, cr.getCount());
        assertTrue(cr.hasNextTransition());

        cr.nextTransition();
        assertEquals(30, cr.getCount());
        assertFalse(cr.hasNextTransition());

        cr.close();
    }

    @Test
    public void testIndices() throws IOException {
        final String basename = FilenameUtils.concat(BASE_TEST_DIR, "102.bin");
        final CountsArchiveWriter writer = new CountsArchiveWriter(basename);
        CountsWriter cw = writer.newCountWriter(0);
        writer.returnWriter(cw);

        cw = writer.newCountWriter(1);

        writer.returnWriter(cw);
        writer.close();
        final CountsArchiveReader reader = new CountsArchiveReader(basename);
        assertEquals(2, reader.getNumberOfIndices());
        assertEquals(2, reader.getIndices().size());
        assertNotNull(reader.getCountReader(0));
        assertNotNull(reader.getCountReader(1));
    }

    @Test
    public void testNames() throws IOException {
        final String basename = FilenameUtils.concat(BASE_TEST_DIR, "103.bin");
        final CountsArchiveWriter writer = new CountsArchiveWriter(basename);
        CountsWriter cw = writer.newCountWriter(0);
        writer.returnWriter(cw);

        cw = writer.newCountWriter(1);

        writer.returnWriter(cw);
        writer.close();
        final CountsArchiveReader reader = new CountsArchiveReader(basename);
        assertEquals(2, reader.getNumberOfIndices());
        assertEquals(2, reader.getIdentifiers().size());
        assertNotNull(reader.getCountReader("0"));
        assertNotNull(reader.getCountReader("1"));
    }
}
