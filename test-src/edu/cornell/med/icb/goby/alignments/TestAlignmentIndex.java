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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.junit.BeforeClass;
import org.junit.AfterClass;

import java.io.IOException;
import java.io.File;

/**
 * @author Fabien Campagne
 *         Date: May 6, 2011
 *         Time: 4:48:24 PM
 */
public class TestAlignmentIndex {
    public void testProblem() throws IOException {
        final AlignmentWriter writer =
                new AlignmentWriter(FilenameUtils.concat(BASE_TEST_DIR, "align-index-error-1"));
        writer.setNumAlignmentEntriesPerChunk(2);

        int position = 1;
        int queryIndex = 0;
        int referenceIndex = 0;

        writer.setAlignmentEntry(queryIndex++, referenceIndex, position++, 30, false, constantQueryLength);
        writer.appendEntry();


        writer.close();

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
