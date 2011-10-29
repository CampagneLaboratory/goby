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

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author campagne
 *         Date: 6/16/11
 *         Time: 4:18 PM
 */
public class TestCachingCountsArchiveReader {
    private static final Log LOG = LogFactory.getLog(TestCountsArchive.class);
    private static final String BASE_TEST_DIR = "test-results/cachingcounts";

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
    public void testCache() throws IOException {


        final String basename = FilenameUtils.concat(BASE_TEST_DIR, "101.bin");
        final CountsArchiveWriter writer = new CountsArchiveWriter(basename);

        CountsWriterI cw = writer.newCountWriter(0, "count-0");
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

        CachingCountsArchiveReader reader = new CachingCountsArchiveReader(basename);
        CountsReaderI countReader = reader.getCountReader("count-0");
        countReader.reposition(10000);
        // read the cached reader again:
        CountsReaderI countReader2 = reader.getCountReader("count-0");
        countReader2.reposition(10000);
    }
}
