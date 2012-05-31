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

package edu.cornell.med.icb.goby.alignments.perms;

import org.apache.commons.io.FileUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author Fabien Campagne
 *         Date: 3/8/12
 *         Time: 3:41 PM
 */
public class TestQueryIndexPermutation {
    @BeforeClass
    public static void setUp() throws IOException {
        FileUtils.forceMkdir(new File("test-results/permutations"));
    }

    @Test
    public void doubleLocal() {

        QueryIndexPermutation perm = new QueryIndexPermutation("test-results/permutations/test-doubleLocal");
        perm.setPruneLimit((byte) 50);
        assertEquals(0, perm.permutate(0));
        assertEquals(0, perm.permutate(0));
        perm.close();
    }

    @Test
    public void simple() {

        QueryIndexPermutation perm = new QueryIndexPermutation("test-results/permutations/test-simple");
        perm.setPruneLimit((byte) 50);
        assertEquals(0, perm.permutate(0));
        assertEquals(1, perm.permutate(1));
        assertEquals(1, perm.permutate(1));
        assertEquals(2, perm.permutate(2));
        assertEquals(3, perm.permutate(3));
        perm.close();

    }

    @Test
    public void multipleChunks() {

        QueryIndexPermutation perm = new QueryIndexPermutation("test-results/permutations/test-multipleChunks-1");
        perm.setPruneLimit((byte) 2);
        assertEquals(0, perm.permutate(0));
        assertEquals(1, perm.permutate(5));
        assertEquals(2, perm.permutate(100));
        assertEquals(1, perm.permutate(5));
        assertEquals(0, perm.permutate(0));
        assertTrue(!perm.isInMap(0));

        assertEquals(2, perm.permutate(100));
        perm.close();
    }


    @Test
    public void multipleChunksVariablePruneLimit() {

        QueryIndexPermutation perm = new QueryIndexPermutation("test-results/permutations/test-multipleChunks-2");
        perm.setPruneLimit((byte) 1);
        assertEquals(0, perm.permutate(0, 3));
        assertEquals(1, perm.permutate(5));
        // 5 has been forgotten and writen to disk already
        assertEquals(2, perm.permutate(100, 2));
        // get -1 5 since it was written to storage already and is requested more than maxOccurences (this indicates a bug in the client).
        assertEquals(-1, perm.permutate(5));
        assertEquals(0, perm.permutate(0, 3));
        // third time's a charm
        assertEquals(0, perm.permutate(0, 3));
        assertEquals(-1, perm.permutate(0, 3));
        assertTrue(!perm.isInMap(0));
        assertTrue(perm.isOnDisk(0));

        assertEquals(2, perm.permutate(100, 2));
        perm.close();
    }
}
