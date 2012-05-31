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
 *         Date: 2/24/12
 *         Time: 4:02 PM
 */
public class TestMedianCombinator   {
    @Test
    public void test1(){
        double[] ps={0.01,0.2,0.4,0.5,0.01};
        MedianCombinator medianC=new MedianCombinator();
        for (double p: ps) {
            medianC.observe(p);
        }
        assertEquals(0.2, medianC.adjust());
    }

    @Test
    public void test2(){
        double[] ps={1,1,1,1,1,0.1};
        MedianCombinator medianC=new MedianCombinator();
        for (double p: ps) {
            medianC.observe(p);
        }
        assertEquals(1.0, medianC.adjust());
    }
}
