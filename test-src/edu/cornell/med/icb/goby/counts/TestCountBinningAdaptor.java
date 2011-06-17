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

import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.*;

/**
 * @author Fabien Campagne
 *         Date: 6/11/11
 *         Time: 12:48 PM
 */
public class TestCountBinningAdaptor {
    CountsReaderTestSupport source;
    private CountBinningAdaptor binAdapter;

    @Before
    public void setUp() throws Exception {
        source = new CountsReaderTestSupport("(10,0)(2,2)(10,0)(4,6)");  //  2*2+4*6/(4+2)=28/6=4.67 is average count over entire range.

    }

    @Test
    public void testWholeReader() throws IOException {
        binAdapter = new CountBinningAdaptor(source, 50);   // bin size overlaps entire source range.
        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(4.67, binAdapter.getAverage(), .1);
        assertEquals(6, binAdapter.getMax(), .1);
        assertFalse(binAdapter.hasNextTransition());

    }

    @Test
    public void testBin5() throws IOException {
        binAdapter = new CountBinningAdaptor(source, 5);   // bin size overlaps a single transition each time.
        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(10, binAdapter.getPosition());
        assertEquals(16, binAdapter.getLength());   // length is 16 because first transition encompasses the first two zon-zero counts.

        assertEquals(4.67, binAdapter.getAverage(), .1);
        assertEquals(6, binAdapter.getMax(), .1);

        assertFalse(binAdapter.hasNextTransition());

    }

    @Test
    public void testBin20() throws IOException {
        binAdapter = new CountBinningAdaptor(source, 20);   // bin size overlaps entire source range.
        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(10, binAdapter.getPosition());
        assertEquals(16, binAdapter.getLength());

        assertEquals(4.67, binAdapter.getAverage(), .1);
        assertEquals(6, binAdapter.getMax(), .1);

        assertFalse(binAdapter.hasNextTransition());

    }


    @Test
    public void testBin1() throws IOException {
        binAdapter = new CountBinningAdaptor(source, 1);   // bin size smaller than transition length
        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(10, binAdapter.getPosition());
        assertEquals(2, binAdapter.getLength());

        assertEquals(2, binAdapter.getAverage(), .1);
        assertEquals(2, binAdapter.getMax(), .1);
        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(22, binAdapter.getPosition());
        assertEquals(4, binAdapter.getLength());

        assertEquals(6, binAdapter.getAverage(), .1);
        assertEquals(6, binAdapter.getMax(), .1);


        assertFalse(binAdapter.hasNextTransition());

    }

  /*  @Test
    public void testBinWithSkipTo() throws IOException {
        String format = "(10,0)(2,2)(10,0)(4,6)";
        CountsReaderTestSupport reader = new CountsReaderTestSupport(format);
        Int2IntMap positionToCounts = new Int2IntArrayMap();
        while (reader.hasNextTransition()) {
            reader.nextTransition();
            for (int i = 0; i < reader.getLength(); i++) {
                positionToCounts.put(reader.getPosition() + i, reader.getCount());
            }

        }
        for (int binSize = 3; binSize < 5; binSize++) {
            for (int skipToPosition = 0; skipToPosition < 30; skipToPosition++) {

                CountsReaderTestSupport readerLocal = new CountsReaderTestSupport(format);
                CountBinningAdaptor bin = new CountBinningAdaptor(readerLocal, binSize);
            //    bin.skipTo(skipToPosition);
                while (bin.hasNextTransition()) {

                    bin.nextTransition();
                    System.out.printf("position=%d count=%g%n",bin.getPosition(), bin.getAverage());
                    double binAverage = 0;
                    int count = 0;
                    for (int j = 0; j < binSize; j++) {
                       final int value = positionToCounts.get(bin.getPosition()+j);
                        binAverage += value;
                        if (value != 0) {
                            count++;
                        }
                    }
                    binAverage /= count;
                    assertEquals(String.format("average differ from expected at skipPosition=%d position=%d binSize=%d %n",
                            skipToPosition, bin.getPosition(), binSize),
                            binAverage, bin.getAverage());
                }
            }
        }
    }       */
}
