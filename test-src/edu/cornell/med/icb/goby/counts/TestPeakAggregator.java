/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.goby.counts;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: May 29, 2009
 *         Time: 1:08:55 PM
 */
public class TestPeakAggregator {
    private static final Log LOG = LogFactory.getLog(TestPeakAggregator.class);
    private static final String BASE_TEST_DIR = "test-results/peaks";

    private final int lengthA = 5;
    private final int lengthB = 100000;
    private final int lengthC = 10;
    private final int lengthD = 1;

    @Test
    public void testPeakFinding() throws IOException {
        final String basename = FilenameUtils.concat(BASE_TEST_DIR, "counts-101.bin");
        writeSomeCounts(basename);
        CountsReader countReader = new CountsReader(new FileInputStream(basename));
        PeakAggregator peakIterator = new PeakAggregator(countReader);
        int count = 0;
        while (peakIterator.hasNext()) {
            final Peak peak = peakIterator.next();
            System.out.println(peak);
            count++;
        }
        assertEquals(2, count);

        //rewind and start again:
        countReader = new CountsReader(new FileInputStream(basename));
        peakIterator = new PeakAggregator(countReader);

        assertTrue(peakIterator.hasNext());
        Peak peak = peakIterator.next();
        assertEquals(5, peak.start);

        assertEquals(lengthB + lengthC, peak.length);
        assertEquals(13, peak.count);

        peak = peakIterator.next();
        assertEquals(lengthA+lengthB+lengthC+lengthA, peak.start);
        assertEquals(lengthD, peak.length);
        assertEquals(10, peak.count);

    }

    private void writeSomeCounts(final String basename) throws IOException {
        final CountsWriter writer = new CountsWriter(new FileOutputStream(basename));
        writer.appendCount(0, lengthA);
        writer.appendCount(1, lengthB);
        writer.appendCount(12, lengthC);
        writer.appendCount(0, lengthA);
        writer.appendCount(10, lengthD);
        writer.appendCount(0, lengthD);
        writer.close();
    }

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
