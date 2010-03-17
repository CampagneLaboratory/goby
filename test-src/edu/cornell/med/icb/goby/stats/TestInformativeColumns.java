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

import org.junit.Test;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

/**
 * Describe class here.
 *
 * @author Kevin Dorff
 */
public class TestInformativeColumns {
    @Test
    public void testInformativeEarly() {
        InformativeColumns ic = new InformativeColumns(3, new InformativeNonZeroNonNaN());
        ic.checkInformative(1.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(3.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(9.0);
        assertTrue(ic.isAllColumnsInformative());
        assertTrue(ic.isColumnInformative(0));
        assertTrue(ic.isColumnInformative(1));
        assertTrue(ic.isColumnInformative(2));
    }

    @Test
    public void testInformativeMid() {
        InformativeColumns ic = new InformativeColumns(3, new InformativeNonZeroNonNaN());
        ic.checkInformative(1.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(0.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(Double.NaN);
        assertFalse(ic.isAllColumnsInformative());

        ic.checkInformative(1.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(0.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(9.0);
        assertFalse(ic.isAllColumnsInformative());

        ic.checkInformative(1.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(3.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(9.0);
        assertTrue(ic.isAllColumnsInformative());

        ic.checkInformative(Double.NaN);
        assertTrue(ic.isAllColumnsInformative());
        ic.checkInformative(Double.NaN);
        assertTrue(ic.isAllColumnsInformative());
        ic.checkInformative(0.0);
        assertTrue(ic.isAllColumnsInformative());

        assertTrue(ic.isColumnInformative(0));
        assertTrue(ic.isColumnInformative(1));
        assertTrue(ic.isColumnInformative(2));
    }

    @Test
    public void testInformativeNever() {
        InformativeColumns ic = new InformativeColumns(3, new InformativeNonZeroNonNaN());
        ic.checkInformative(1.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(0.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(Double.NaN);
        assertFalse(ic.isAllColumnsInformative());

        ic.checkInformative(1.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(0.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(9.0);
        assertFalse(ic.isAllColumnsInformative());

        ic.checkInformative(1.0);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(Double.NaN);
        assertFalse(ic.isAllColumnsInformative());
        ic.checkInformative(9.0);
        assertFalse(ic.isAllColumnsInformative());

        assertTrue(ic.isColumnInformative(0));
        assertFalse(ic.isColumnInformative(1));
        assertTrue(ic.isColumnInformative(2));
    }
}
