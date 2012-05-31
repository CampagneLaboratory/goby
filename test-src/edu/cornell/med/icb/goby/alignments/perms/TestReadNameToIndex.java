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

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 3/10/12
 *         Time: 5:02 PM
 */
public class TestReadNameToIndex {
    @Test
    public void testSimple(){
        ReadNameToIndex rnti=new ReadNameToIndex("test-results/permutations/rnti-simple-1.tsv");
        assertEquals(0, rnti.getQueryIndex("name-0", 2));
        assertEquals(0, rnti.getQueryIndex("name-0", 2));
        assertEquals(1, rnti.getQueryIndex("name-0", 2));
        assertEquals(1, rnti.getQueryIndex("name-0", 2));
        assertEquals(2, rnti.getQueryIndex("name-0", 2));
        assertEquals(2, rnti.getQueryIndex("name-0", 2));
        assertEquals(3, rnti.getQueryIndex("name-0", 2));
        assertEquals(3, rnti.getQueryIndex("name-0", 2));
    }

    @Test
    public void testVarMaxOcc(){
        ReadNameToIndex rnti=new ReadNameToIndex("test-results/permutations/rnti-simple-1.tsv");
        assertEquals(0, rnti.getQueryIndex("name-0", 0));
        assertEquals(1, rnti.getQueryIndex("name-0", 0));
        assertEquals(2, rnti.getQueryIndex("name-1", 1));
        assertEquals(3, rnti.getQueryIndex("name-1", 1));
        assertEquals(4, rnti.getQueryIndex("name-2", 2));
        assertEquals(4, rnti.getQueryIndex("name-2", 2));
        assertEquals(5, rnti.getQueryIndex("name-3", 3));
        assertEquals(5, rnti.getQueryIndex("name-3", 3));
        assertEquals(5, rnti.getQueryIndex("name-3", 3));
        assertEquals(6, rnti.getQueryIndex("name-3", 3));
    }
}
