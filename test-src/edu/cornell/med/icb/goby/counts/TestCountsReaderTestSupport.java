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

import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: May 21, 2011
 *         Time: 11:15:57 AM
 */
public class TestCountsReaderTestSupport {
    private CountsReaderTestSupport countsReader;
    private int[] lengths;
    private int[] counts;

    @Test
    public void iterate() throws IOException {

        while (countsReader.hasNextTransition()) {
            countsReader.nextTransition();

        }
    }

    @Test
    public void validate() throws IOException {

        assertTrue(countsReader.hasNextTransition());
        assertEquals(0, countsReader.getPosition());

        countsReader.nextTransition();

        assertEquals(5, countsReader.getPosition());
        assertEquals(0, countsReader.getCount());
        assertEquals(5, countsReader.getLength());

        countsReader.nextTransition();
        assertEquals(6, countsReader.getPosition());
        assertEquals(1, countsReader.getCount());
        assertEquals(1, countsReader.getLength());

        countsReader.nextTransition();
        assertEquals(26, countsReader.getPosition());
        assertEquals(0, countsReader.getCount());
        assertEquals(20, countsReader.getLength());

        countsReader.nextTransition();
        assertEquals(56, countsReader.getPosition());
        assertEquals(20, countsReader.getCount());
        assertEquals(30, countsReader.getLength());

        countsReader.nextTransition();
        assertEquals(1056, countsReader.getPosition());
        assertEquals(0, countsReader.getCount());
        assertEquals(1000, countsReader.getLength());

        assertFalse(countsReader.hasNextTransition());

    }


    @Before
    public void setup() {
        counts = new int[]{0, 1, 0, 20, 0};
        lengths = new int[]{5, 1, 20, 30, 1000};
        countsReader = new CountsReaderTestSupport(lengths, counts);
    }
}
