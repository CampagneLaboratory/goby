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

package edu.cornell.med.icb.goby.modes;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
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
 *         Date: 4/29/12
 *         Time: 11:06 AM
 */
public class TestSortMode {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestSortMode.class);
    private static final String BASE_TEST_DIR = "test-results/sort-alignment";

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
      //  FileUtils.forceDeleteOnExit(new File(BASE_TEST_DIR));
    }

    @Test
    // check that large-sort can sort a small alignment (only one split):
    public void sortSmall() throws IOException {
        SortMode mode = new SortMode();
        mode.setInput("test-data/alignment-hybrid-codec/EJOYQAZ-small.header");
        mode.setOutput(FilenameUtils.concat(BASE_TEST_DIR, "EJOYQAZ-small-sorted"));
        mode.setNumThreads(1);
        mode.setSplitSize(1 * 1024 * 1024);
        mode.execute();
        // quick way to check the entries and header files were created. Since it is a checksum, if the test fails, it is possible
        // content is correct, but written sligthly differently.
        assertEquals(1947630632, FileUtils.checksumCRC32(new File(FilenameUtils.concat(BASE_TEST_DIR, "EJOYQAZ-small-sorted.entries")))) ;
        assertEquals(2375098082L, FileUtils.checksumCRC32(new File(FilenameUtils.concat(BASE_TEST_DIR, "EJOYQAZ-small-sorted.header")))) ;
    }
}
