/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.stats;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

/**
 * Test InformativeNonZeroNonNaN.
 *
 * @author Kevin Dorff
 */
public class TestInformativeNonZeroNonNaN {
    @Test
    public void testInformative() {
        final InformativeDouble informativeObj = new InformativeNonZeroNonNaN();
        assertFalse(informativeObj.isInformative(Double.NaN));
        assertFalse(informativeObj.isInformative(0));
        assertTrue(informativeObj.isInformative(-10));
        assertTrue(informativeObj.isInformative(10));
        assertTrue(informativeObj.isInformative(Double.MAX_VALUE));
        assertTrue(informativeObj.isInformative(Double.MIN_VALUE));
        assertTrue(informativeObj.isInformative(Double.POSITIVE_INFINITY));
        assertTrue(informativeObj.isInformative(Double.NEGATIVE_INFINITY));
    }
}
