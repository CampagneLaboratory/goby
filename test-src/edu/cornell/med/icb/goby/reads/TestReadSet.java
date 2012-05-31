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

package edu.cornell.med.icb.goby.reads;

import org.apache.commons.io.FileUtils;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Jun 3, 2009
 *         Time: 3:55:50 PM
 */
public class TestReadSet {
    @BeforeClass
    public static void initializeTestDirectory() throws IOException {
        FileUtils.forceMkdir(new File("test-results/read-sets"));
    }

    @Test
    public void testSave() throws IOException {
        final ReadSet set1 = new ReadSet();
        set1.smallestStoredMultiplicity(0);
        set1.add(2);
        set1.add(0);
        set1.add(1);

        set1.add(3);
        set1.save("test-results/read-sets/set", "1");

        final ReadSet set2 = new ReadSet();
        set2.load("test-results/read-sets/set", "1");
        assertEquals(set1.size(), set2.size());
        assertTrue(set2.contains(0));
        assertTrue(set2.contains(1));
        assertTrue(set2.contains(2));
        assertTrue(set2.contains(3));
        assertFalse(set2.contains(4));
    }

    @Test
    public void testMultiplicity() throws IOException {
        final ReadSet set1 = new ReadSet();
        set1.add(2, 1001);
        set1.add(0, 1);
        set1.add(1, 5);

        set1.add(3, 13);
        set1.save("test-results/read-sets/set", "2");

        final ReadSet set2 = new ReadSet();
        set2.load("test-results/read-sets/set", "2");
        assertEquals(set1.size(), set2.size());
        assertTrue(set2.contains(0));
        assertTrue(set2.contains(1));
        assertTrue(set2.contains(2));
        assertTrue(set2.contains(3));
        assertFalse(set2.contains(4));

        assertEquals(1001, set2.getMultiplicity(2));
        assertEquals(1, set2.getMultiplicity(0));
        assertEquals(5, set2.getMultiplicity(1));
        assertEquals(13, set2.getMultiplicity(3));
    }

    @Test
    public void testSmallestMultiplicity() throws IOException {
        final ReadSet set1 = new ReadSet();
        set1.smallestStoredMultiplicity(1);
        set1.add(2, 1001);
        set1.add(0, 1);
        set1.add(1, 5);

        set1.add(3, 13);
        set1.save("test-results/read-sets/set", "2");

        final ReadSet set2 = new ReadSet();
        set2.load("test-results/read-sets/set", "2");
        assertEquals(set1.size(), set2.size());
        assertTrue(set2.contains(0));
        assertTrue(set2.contains(1));
        assertTrue(set2.contains(2));
        assertTrue(set2.contains(3));
        assertFalse(set2.contains(4));

        assertEquals(1001, set2.getMultiplicity(2));
        assertEquals(1, set2.getMultiplicity(0));
        assertEquals(5, set2.getMultiplicity(1));
        assertEquals(13, set2.getMultiplicity(3));
        assertEquals(0, set2.getMultiplicity(50));
    }
}
