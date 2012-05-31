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
 *         Date: 2/23/12
 *         Time: 2:20 PM
 */
public class TestQFast {

    @Test
    public void test1() {
        double product=1;
        product*=0.01;
        product*=0.01;
        assertEquals(0.0001, QFast.qfast(2, product), 1e-6);

        product=1;
        product*=0.1;
        product*=0.001;
        assertEquals(0.0001, QFast.qfast(2, product), 1e-6);

        product=1;
        product*=0.1;
        product*=0.001;
        product*=0.5;
        product*=0.5;
        //System.out.printf("product="+product);
        assertEquals(0.001693524214160628, QFast.qfast(4, product), 1e-6);
    }


}
