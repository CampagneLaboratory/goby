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

import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.CumulativeArrayKeeperFast;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 *
 * @Author nyasha
 * Date: 3/2/12
 * Time: 3:59 PM
 */

public class TestCumulativeArrayKeeper {

    @Test
    public void testCase1() {

        CumulativeArrayKeeperFast cumCtest= new CumulativeArrayKeeperFast(5);
        cumCtest.addToRight(5);
        assertEquals("[\t5\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(10);
        cumCtest.addToRight(2);
        cumCtest.addToRight(8);
        cumCtest.addToRight(10);
        assertEquals("[\t5\t15\t17\t25\t35\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(6);
        assertEquals("[\t10\t12\t20\t30\t36\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(3);
        assertEquals("[\t2\t10\t20\t26\t29\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(8);
        assertEquals("[\t8\t18\t24\t27\t35\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(5);
        assertEquals("[\t10\t16\t19\t27\t32\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(9);
        assertEquals("[\t6\t9\t17\t22\t31\t]", cumCtest.getCumulativeArray());

    }

    @Test
    public void testCase2(){
        CumulativeArrayKeeperFast cumCtest= new CumulativeArrayKeeperFast(3);
        cumCtest.addToRight(0);
        assertEquals("[\t0\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(5);
        assertEquals("[\t0\t5\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(0);
        assertEquals("[\t0\t5\t5\t]", cumCtest.getCumulativeArray());
    }

    @Test
    public void testCase3(){
        CumulativeArrayKeeperFast cumCtest= new CumulativeArrayKeeperFast(3);
        cumCtest.addToRight(0);
        assertEquals("[\t0\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(0);
        assertEquals("[\t0\t0\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(0);
        assertEquals("[\t0\t0\t0\t]", cumCtest.getCumulativeArray());
    }

    @Test
    public void testCase4(){
        CumulativeArrayKeeperFast cumCtest= new CumulativeArrayKeeperFast(1);
        cumCtest.addToRight(0);
        assertEquals("[\t0\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(5);
        assertEquals("[\t5\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(456);
        assertEquals("[\t456\t]", cumCtest.getCumulativeArray());
        cumCtest.addToRight(0);
        assertEquals("[\t0\t]", cumCtest.getCumulativeArray());
    }


}
