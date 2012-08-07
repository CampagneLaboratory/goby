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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import junit.framework.TestCase;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: 7/27/12
 *         Time: 12:06 PM
 */
public class RandomAccessAnnotationsTest extends TestCase {

    public void testOverlap() {
        RandomAccessAnnotations finder = new RandomAccessAnnotations();
        finder.addAnnotation("element-1", "chrX", 10, 310);          // gene starts at 10, ends at 310
        finder.addAnnotation("element-1", "chrX", 1000, 1001);          // gene starts at 1000, ends at 1001
        assertNotNull(finder.find("chrX", 11, 200));
        assertNotNull(finder.find("chrX", 10, 310));
        assertNull(finder.find("chrX", 8, 9));
        assertNull(finder.find("chrX", 311, 312));
        assertNotNull(finder.find("chrX", 1000, 1001));
        assertNull(finder.find("chrX", 1002, 1003));
    }

    public void testLoad() throws IOException {
        RandomAccessAnnotations finder = new RandomAccessAnnotations();
        finder.loadAnnotations("/Users/campagne/IdeaProjects/gobyweb/exon-annotations.tsv");
        assertNotNull(finder.find("1", 11, 200));
        assertNotNull(finder.find("1", 10, 310));
        assertNull(finder.find("1", 8, 9));
        assertNull(finder.find("1", 311, 312));
        assertNotNull(finder.find("1", 1000, 1001));
        assertNull(finder.find("1", 1002, 1003));
    }
}
