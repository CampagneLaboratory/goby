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

package edu.cornell.med.icb.goby.modes;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * Make sure read-quality-stats is returning stats if the file has quality scores.
 */
public class TestReadQualityStatsMode {
    @Test
    public void testHasStats() throws IOException {
        final ReadQualityStatsMode readQuality = new ReadQualityStatsMode();
        readQuality.setSampleFraction(1.0d);
        readQuality.addInputFile(new File("test-data/compact-reads/five-with-quality.compact-reads"));
        final File outputFile = new File("test-results/five-with-quality.compact-reads.read-qual-stats.tsv");
        readQuality.setOutputFile(outputFile);
        readQuality.execute();
        assertEquals("All reads should have been observed", 5, readQuality.getNumberOfObservedReads());
        assertEquals("No reads should have been skipped", 0, readQuality.getNumberOfSkippedReads());
        final List lines = FileUtils.readLines(outputFile);
        assertEquals("Wrong number of stats lines in output file", 39, lines.size());
    }

    @Test
    public void testHasStatsWithSampling() throws IOException {
        final ReadQualityStatsMode readQuality = new ReadQualityStatsMode();
        readQuality.setSampleFraction(0.5d);  // Keep about half of the reads
        readQuality.addInputFile(new File("test-data/compact-reads/five-with-quality.compact-reads"));
        final File outputFile = new File("test-results/five-with-quality.compact-reads.read-qual-stats.tsv");
        readQuality.setOutputFile(outputFile);
        readQuality.execute();
        assertTrue("All reads should have been observed", readQuality.getNumberOfObservedReads() < 5);
        assertTrue("No reads should have been skipped", readQuality.getNumberOfSkippedReads() > 0);
        final List lines = FileUtils.readLines(outputFile);
        assertEquals("Wrong number of stats lines in output file", 39, lines.size());
    }
}
