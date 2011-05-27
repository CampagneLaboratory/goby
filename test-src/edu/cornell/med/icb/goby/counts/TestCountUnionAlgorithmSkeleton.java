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

import static junit.framework.Assert.assertEquals;

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
        CountsReaderI reader0 = new CountsReaderTestSupport("(1,0)(3,3)(2,-3)");
        CountsReaderI reader1 = new CountsReaderTestSupport("(6,1)");

        String expected ="(1,1)(3,3)(2,-3)";

        UnionAlgorithmSkeleton unionAlgo=new UnionAlgorithmSkeleton(reader0, reader1);

        CountsReaderI reader = unionAlgo;
        MutableString result=new MutableString();

        while (reader.hasNextTransition()) {
            reader.nextTransition();
            result.append(String.format("(%d,%d)", reader.getPosition(),
                    reader.getCount()));

        }
        System.out.println("result: "+result);
        assertEquals("the result of the union must match", expected,result.toString());
    }
}
