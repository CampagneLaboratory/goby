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

package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 3/1/12
 *         Time: 11:34 AM
 */
public class TestSitesInFixedWindow {
    @Test
    public void testAdd() {
        SitesInFixedWindow window = new SitesInFixedWindow(5);
        window.add(0, 0);
        assertEquals(1, window.n());
        window.add(0, 1);
        assertEquals(2, window.n());
        window.add(0, 2);
        assertEquals(3, window.n());
        window.add(0, 3);
        window.add(0, 4);
        window.add(0, 5);
        assertEquals(6, window.n());
        window.add(0, 6);
        assertEquals(6, window.n());
        window.add(0, 7);
        assertEquals(6, window.n());
    }

    @Test
    public void testPrune() {
        SitesInFixedWindow window = new SitesInFixedWindow(5);
        window.add(0, 0);

        window.add(0, 20);
        assertEquals(1, window.n());

        window.add(0, 21);
        window.add(0, 22);
        window.add(0, 23);
        assertEquals(4, window.n());

    }
}
