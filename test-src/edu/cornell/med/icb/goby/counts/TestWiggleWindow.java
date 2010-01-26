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

package edu.cornell.med.icb.goby.counts;

import org.junit.Test;

import java.io.PrintWriter;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import static junit.framework.Assert.assertEquals;

/**
 * Describe class here.
 *
 * @author Kevin Dorff
 */
public class TestWiggleWindow {

    @Test
    public void testWiggleWindow() throws IOException {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        final WiggleWindow wiggleWindow = new WiggleWindow(new PrintWriter(baos), 10, 100);
        wiggleWindow.addData(0, 3, 3);
        wiggleWindow.addData(4, 3, 4);
        wiggleWindow.addData(11, 3, 5);
        wiggleWindow.addData(16, 6, 6);
        wiggleWindow.addData(80, 34, 7);
        wiggleWindow.addData(130, 24, 8);
        wiggleWindow.finish();
        baos.close();
        String actual = baos.toString();
        String expected = String.format("1 3%n12 5%n22 6%n81 7%n91 7%n");
        assertEquals(expected, actual);
    }

}
