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

package edu.cornell.med.icb.goby.algorithmic.compression;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 1/15/12
 *         Time: 11:15 AM
 */
public class TestMoveToFrontCoder {
    @Test
    public void testEncode() {
        MoveToFrontCoder coder = new MoveToFrontCoder(4);

        assertEquals(0, coder.encode(0));
        assertEquals(0, coder.encode(0));
        assertEquals(1, coder.encode(1));
        assertEquals(0, coder.encode(1));
        assertEquals(3, coder.encode(3));
    }
}
