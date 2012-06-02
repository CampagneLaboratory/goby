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

import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * Test of sample quality encodings.
 */
public class TestSampleQualityScoresMode {
    @Test
    public void testFasta() throws IOException {
        final SampleQualityScoresMode sqs = new SampleQualityScoresMode();
        sqs.addInputFilename("test-data/sample-qual-scores/30reads.fa");
        sqs.execute();
        final List<String> encodings = sqs.getLikelyEncodings();
        assertEquals("wrong number of encodings", 1, encodings.size());
        assertEquals("wrong encoding", "fasta", encodings.get(0));
    }

    @Test
    public void testFastq() throws IOException {
        final SampleQualityScoresMode sqs = new SampleQualityScoresMode();
        sqs.addInputFilename("test-data/sample-qual-scores/30reads.fq");
        sqs.execute();
        final List<String> encodings = sqs.getLikelyEncodings();
        assertEquals("wrong number of encodings", 1, encodings.size());
        assertEquals("wrong encoding", "Illumina/Solexa", encodings.get(0));
    }
}
