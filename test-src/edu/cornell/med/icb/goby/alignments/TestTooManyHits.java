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

package edu.cornell.med.icb.goby.alignments;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: May 6, 2009
 *         Time: 1:46:14 PM
 */
public class TestTooManyHits {
    @Test
    public void testTooManyHits() throws IOException {
        final AlignmentTooManyHitsWriter tmhWriter =
                new AlignmentTooManyHitsWriter("test-results/alignments/align-101-tmh", 4);

        // tmhWriter will only write entries if numHits > thresh
        tmhWriter.getNewAmbiguousLocation().setAtLeastNumberOfHits(5);
        tmhWriter.getNewAmbiguousLocation().setQueryIndex(0);
        tmhWriter.append();
        tmhWriter.getNewAmbiguousLocation().setAtLeastNumberOfHits(15);
        tmhWriter.getNewAmbiguousLocation().setQueryIndex(12);
        tmhWriter.append();
        tmhWriter.close();

        final AlignmentTooManyHitsReader tmhReader =
                new AlignmentTooManyHitsReader("test-results/alignments/align-101-tmh");
        assertTrue("query sequence 0 must be found", tmhReader.isQueryAmbiguous(0));
        /*
        There are three ambiguity-related values involved in the following test.
         A- tmhWriter.alignerThreshold
         B- query sequence 0's numHits value
         C- argument #2 supplied to isQueryAmbiguous

         If A >= B, then query sequence 0 is not written.  The test should be assertFalse
         If A < C < B, then query sequence 0 is written, and test should be assertTrue
         If A < B <= C, then query sequence 0 is written, but *counter-intuitively* isQueryAmbiguous returns true
         */
        assertTrue("query sequence 0 must be found at k=3", tmhReader.isQueryAmbiguous(0, 3));
        assertTrue("query sequence 0 *counter-intuitively* must be found at k=20 because 20 is larger than the aligner threshold", tmhReader.isQueryAmbiguous(0, 20));
        assertFalse("query sequence 1 must NOT be found", tmhReader.isQueryAmbiguous(1));
        assertTrue("query sequence 12 must be found", tmhReader.isQueryAmbiguous(12));
        assertFalse("query sequence 13 must NOT be found", tmhReader.isQueryAmbiguous(13));
        assertFalse("query sequence 1100239028 must NOT be found", tmhReader.isQueryAmbiguous(1100239028));

    }
}
