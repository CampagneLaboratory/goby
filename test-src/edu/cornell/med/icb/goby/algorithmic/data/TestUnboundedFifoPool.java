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

package edu.cornell.med.icb.goby.algorithmic.data;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

/**
 * @author Fabien Campagne
 *         Date: Apr 30, 2011
 *         Time: 1:36:00 PM
 */
public class TestUnboundedFifoPool {
    @Test
    public void simple() {
        UnboundedFifoPool<Integer> tester = new UnboundedFifoPool<Integer>();
        assertTrue(tester.isEmpty());
        tester.add(1);
        assertEquals((Integer) 1, tester.remove());
        tester.add(2);
        tester.add(3);

        assertEquals((Integer) 2, tester.remove());
        assertEquals((Integer) 3, tester.remove());
        assertTrue(tester.isEmpty());
        tester.add(4);
        tester.add(5);
        tester.add(6);
        assertEquals((Integer) 4, tester.remove());
        assertEquals((Integer) 5, tester.remove());
        assertEquals((Integer) 6, tester.remove());
        assertTrue(tester.isEmpty());
    }


}
