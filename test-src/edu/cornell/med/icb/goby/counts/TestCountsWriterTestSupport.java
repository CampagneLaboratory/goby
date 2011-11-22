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

import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;


/**
 * @author Fabien Campagne
 *         Date: 10/29/11
 *         Time: 1:57 PM
 */
public class TestCountsWriterTestSupport {
    @Test
    public void simple() throws IOException {
        CountsWriterTestSupport writer = new CountsWriterTestSupport();
        writer.appendCount(10, 4);
        writer.appendCount(1, 0);
        writer.appendCount(2, 8);
        assertEquals("initial-count=0 (c=10,l=4)(c=1,l=0)(c=2,l=8)", writer.countsAsText());
    }

    @Test
    public void repeatingCounts() throws IOException {
        final CountsWriterTestSupport writer = new CountsWriterTestSupport();
        System.out.println("Hello");
        writer.appendCount(10, 4);
        try {
            writer.appendCount(10, 1);
            System.out.println("Failed");
            throw new InternalError("assertion must have been raised");

        } catch (AssertionError e) {
            // OK

        }

    }
}
