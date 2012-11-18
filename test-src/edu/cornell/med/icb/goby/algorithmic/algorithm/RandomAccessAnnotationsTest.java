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

import edu.cornell.med.icb.goby.algorithmic.data.Interval;
import it.unimi.dsi.lang.MutableString;
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
        finder.loadAnnotations("data/biomart-mouse-exons-ensembl57-genes-MM9.txt");
        /// chr17	-1	ENSMUSG00000077184	ENSMUSE00000664149	24,160,167	24,160,293
///      chr4	1	ENSMUSG00000083154	ENSMUSE00000717873	    24,165,247	24,165,247
        ///   chr10	1	ENSMUSG00000086148	ENSMUSE00000794298	24,167,296	24,172,463

        final Interval interval1 = finder.find("chr4", 24165247, 24165247);
        assertNotNull(interval1);
        final int chr4Index = finder.references.getInt(new MutableString("chr4"));

        assertEquals(chr4Index, interval1.referenceIndex);

        final Interval interval2 = finder.find("chr11",  11838488, 11843885);
        assertNotNull(interval2);
        final int chr11Index = finder.references.getInt(new MutableString("chr11"));

        assertEquals(chr11Index, interval2.referenceIndex);

    }
}
