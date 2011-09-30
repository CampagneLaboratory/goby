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

package edu.cornell.med.icb.goby.util.motifs;

import it.unimi.dsi.util.CircularCharArrayBuffer;
import org.junit.Test;

import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertTrue;

/**
 * @author campagne
 *         Date: 9/29/11
 *         Time: 12:32 PM
 */
public class TestSubSequenceMotifMatcher {
    @Test
    public void testMotifLength1() {

        SubSequenceMotifMatcher matcher = new SubSequenceMotifMatcher("C");
        matcher.accept('A');
        assertFalse(matcher.matched());
        matcher.accept('C');
        assertTrue(matcher.matched());
        matcher.accept('T');
        assertFalse(matcher.matched());
        matcher.accept('G');
        assertFalse(matcher.matched());
        matcher.accept('F');
        assertFalse(matcher.matched());
        matcher.accept('C');
        assertTrue(matcher.matched());
        matcher.newSequence();
        matcher.accept('C');
        assertTrue(matcher.matched());

    }

    @Test
    public void testMotifLength2() {

        SubSequenceMotifMatcher matcher = new SubSequenceMotifMatcher("CG");
        matcher.accept('A');
        assertFalse(matcher.matched());
        matcher.accept('C');
        assertFalse(matcher.matched());
        matcher.accept('T');
        assertFalse(matcher.matched());
        matcher.accept('C');
        assertFalse(matcher.matched());
        matcher.accept('G');
        assertTrue(matcher.matched());
        matcher.accept('G');
        assertFalse(matcher.matched());

    }
}
