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

package edu.cornell.med.icb.goby.counts;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntListIterator;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 10/29/11
 *         Time: 2:33 PM
 */
public class TestCountsWriterHelper {

    @Test
    public void testCase1() throws IOException {
        CountsWriterTestSupport delegate = new CountsWriterTestSupport();
        CountWriterHelper helper = new CountWriterHelper(delegate);

        String input1 = "(1,4)(2,4)(3,4)(4,1)(5,1)(6,1)(9,3)(10,1)";
        String expectedOutput = "initial-count=0 (c=0,l=1)(c=4,l=3)(c=1,l=3)(c=0,l=2)(c=3,l=1)(c=1,l=1)(c=0,l=1)";
        exerciseDelegate(helper, input1);
        assertEquals("counts must match expected, format (length,count):", expectedOutput, delegate.countsAsText());
    }
  //   @Test
    public void testCase2() throws IOException {
        CountsWriterTestSupport delegate = new CountsWriterTestSupport(4);
        CountWriterHelper helper = new CountWriterHelper(delegate);

        String input1 = "(0,4)(1,4)(2,4)(8,3)(11,12)(14,3)(15,2)(16,0)";
        String expectedOutput = "initial-count=0 (c=4,l=3)(c=0,l=5)(c=3,l=1)(c=0,l=2)(c=2,l=1)(c=0,l=2)(c=3,l=1)(c=2,l=1)(c=0,l=1)";
        exerciseDelegate(helper, input1);
        assertEquals("counts must match expected, format (length,count):", expectedOutput, delegate.countsAsText());
    }

    private void exerciseDelegate(CountWriterHelper helper, String input1) throws IOException {
        IntList positions = new IntArrayList();
        IntList counts = new IntArrayList();
        parse(input1, positions, counts);
        IntListIterator positionIt = positions.iterator();
        for (IntListIterator iterator = counts.iterator(); iterator.hasNext(); ) {

            int count = iterator.nextInt();
             positionIt.hasNext();
            int position = positionIt.nextInt();
            helper.appendCountAtPosition(count, position);
        }
        helper.close();
    }

    public int parse(String format, IntList positions, IntList counts) {
        final String[] tokens = format.split("[() ]+");

        for (int i = 0; i < tokens.length; i++) {
            String token = tokens[i];
            if (token.length() > 0) {
                String[] t = token.split(",");
                positions.add(Integer.parseInt(t[0]));
                counts.add(Integer.parseInt(t[1]));
            }
        }
        return positions.size();

    }

}
