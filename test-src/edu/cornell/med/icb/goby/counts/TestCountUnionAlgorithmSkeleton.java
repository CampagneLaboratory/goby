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

import it.unimi.dsi.lang.MutableString;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;


/**
 * Skeleton for testing union algorithm.
 *
 * @author Fabien Campagne
 *         Date: 5/26/11
 *         Time: 9:47 PM
 */
public class TestCountUnionAlgorithmSkeleton {
    @Test
    /**
     * https://docs.google.com/document/d/1VtHTd3Bq9iFC-pAISVtnEBt8EbkWZ9Ot356pjRE5poc/edit?hl=en_US#
     */
    public void testCase1() throws IOException {


        assertEquals("the result of the union must match", "(0,1)(1,4)(4,1)(6,0)", unionPositionCount("(1,0)(3,3)(2,0)", "(6,1)"));
    }


    @Test
    /**
     * https://docs.google.com/document/d/1VtHTd3Bq9iFC-pAISVtnEBt8EbkWZ9Ot356pjRE5poc/edit?hl=en_US#
     */
    public void testCase2() throws IOException {

       assertEquals("the result of the union must match", "(0,0)(1,1)(2,0)(4,1)(6,2)(8,1)(9,2)(10,1)(12,2)(13,1)(14,0)(16,1)(17,0)",
                unionPositionCount("(1,0)(1,1)(4,0)(2,1)(1,0)(1,1)(2,0)(1,1)(3,0)(1,1)",
                        "(4,0)(10,1)"));

        assertEquals("the result of the union must match", "(1,0)(1,1)(2,0)(2,1)(2,2)(1,1)(1,2)(2,1)(1,2)(1,1)(2,0)(1,1)(0,0)",
                unionLengthCount("(1,0)(1,1)(4,0)(2,1)(1,0)(1,1)(2,0)(1,1)(3,0)(1,1)",
                        "(4,0)(10,1)"));
    }

    @Test
    /**
     * https://docs.google.com/document/d/1VtHTd3Bq9iFC-pAISVtnEBt8EbkWZ9Ot356pjRE5poc/edit?hl=en_US#
     */
    public void testCase3() throws IOException {

   //     assertEquals("(1,0)(3,1)", unionLengthCount("(1,0)(3,1)"));
        assertEquals("(0,0)(1,1)(4,0)", unionPositionCount("(1,0)(3,1)"));

    }


    public static String unionLengthCount(final String... formats) throws IOException {
        CountsReaderI readers[] = new CountsReaderI[formats.length];
        int i = 0;
        for (String format : formats) {
            readers[i] = new CountsReaderTestSupport(format);
            i++;
        }
        CountsAggregatorI orIterator;
        MutableString result = new MutableString();
        orIterator = new UnionAlgorithmSkeleton(readers);
        while (orIterator.hasNextTransition()) {
            orIterator.nextTransition();

            result.append(String.format("(%d,%d)", orIterator.getLength(), orIterator.getCount()));

        }
        result.append(String.format("(%d,%d)", 0, orIterator.getCount()));
        return result.toString();
    }

    public static   String unionPositionCount(final String... formats) throws IOException {
        CountsReaderI readers[] = new CountsReaderI[formats.length];
        int i = 0;
        for (String format : formats) {
            readers[i] = new CountsReaderTestSupport(format);
            i++;
        }
        CountsAggregatorI orIterator;
        MutableString result = new MutableString();
        orIterator = new UnionAlgorithmSkeleton(readers);
        while (orIterator.hasNextTransition()) {
            orIterator.nextTransition();

            result.append(String.format("(%d,%d)", orIterator.getPosition(), orIterator.getCount()));
            System.out.println(result);
        }
        // last position must be correctly updated after hasNext returns false:
        result.append(String.format("(%d,%d)", orIterator.getPosition(), 0));
        return result.toString();
    }


}
