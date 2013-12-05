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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 2/20/12
 *         Time: 4:44 PM
 */
public class TestFenwickTree {
    @Test
    public void testIncrement() {
        FenwickTree tree = new FenwickTree(4);
        tree.incrementCount(0);
        assertEquals(1, tree.getCumulativeCount(0));
        tree.incrementCount(0);
        assertEquals(2, tree.getCumulativeCount(0));
        tree.incrementCount(0);
        assertEquals(3, tree.getCumulativeCount(0));
        tree.incrementCount(1);
        assertEquals(3, tree.getCumulativeCount(0));
        assertEquals(4, tree.getCumulativeCount(1));
        tree.incrementCount(3);
        assertEquals(3, tree.getCumulativeCount(0));
        assertEquals(4, tree.getCumulativeCount(1));
        assertEquals(4, tree.getCumulativeCount(2));
        assertEquals(5, tree.getCumulativeCount(3));
    }
<<<<<<< HEAD


=======
    /*
    @Test
    public void testFind() {
        FenwickTree tree = new FenwickTree(4);
        add(tree, 4, 0);  //4 (cumulativeCount up to, including index)
        add(tree, 2, 1);  //6
        add(tree, 1, 2);  //7
        add(tree, 6, 3);  //13
        FenwickTree.answer result = tree.createAnswer();
        tree.find(6, result);
        assertEquals(6, result.cumulativeCount);
        assertEquals(1, result.symbolIndex);
        assertEquals(0, tree.find(4, result));
        assertEquals(2, tree.find(7, result));
        assertEquals(3, tree.find(13, result));
        assertEquals(-1, tree.find(12, result));

        tree.find(6, result);
        assertEquals(6, result.cumulativeCount);
        assertEquals(7, result.nextLargerCumulativeCount);

        tree.find(7, result);
               assertEquals(7, result.cumulativeCount);
               assertEquals(13, result.nextLargerCumulativeCount);



    }
       */
>>>>>>> 413d54e52a172d34cf0a704e1d602d291074d7d9
    private void add(FenwickTree tree, int num, int index) {
        for (int i = 0; i < num; i++) {
            tree.incrementCount(index);
        }
    }


}
