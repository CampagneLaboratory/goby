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

package edu.cornell.med.icb.goby.util;

import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import static org.junit.Assert.assertEquals;
/**
 * @author Fabien Campagne
 *         Date: Apr 1, 2011
 *         Time: 1:07:09 AM
 */
public class TestFiles {
    public void assertEquals(File file1, File file2) throws FileNotFoundException {
        MutableString string1 = new MutableString();
        MutableString string2 = new MutableString();
        string1 = readFile1(file1);
        string2 = readFile1(file2);
        org.junit.Assert.assertEquals(string1.toString(), string2.toString());
    }
    public void assertEquals(int v1, int v2) {
              org.junit.Assert.assertEquals(v1,v2);
        }

    private MutableString readFile1(File file1) throws FileNotFoundException {
        LineIterator it = new LineIterator(new FastBufferedReader(new FileReader(file1)));
        MutableString result = new MutableString();
        while (it.hasNext()) {
            MutableString mutableString = it.next();
            result.append(mutableString);
            result.append("\n");

        }
        return result;
    }
}
