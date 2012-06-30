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

package edu.cornell.med.icb.goby.algorithmic.compression;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.io.FastByteArrayOutputStream;
import it.unimi.dsi.io.DebugInputBitStream;
import it.unimi.dsi.io.DebugOutputBitStream;
import it.unimi.dsi.io.InputBitStream;
import it.unimi.dsi.io.OutputBitStream;
import junit.framework.TestCase;
import org.junit.Test;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: 6/30/12
 *         Time: 3:02 PM
 */
public class TestFastArithmeticCoderOrder1 extends TestCase {
  //  final int[] list1 = {0,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,1, 2, 3, 4, 3, 1, 2, 1, 2,3, 1, 2, 1, 2, 1};
   final int[] list1 = {0,1,2,3,4};
   // final int[] list1 = {0,1,0,12};

    @Test
    public void testOrder0() throws IOException {

        IntArrayList list = new IntArrayList(list1);
        int numSymbols = 5;
        FastArithmeticCoderI order1 = new FastArithmeticCoder(numSymbols);
        FastArithmeticDecoderI decoder = new FastArithmeticDecoder(numSymbols);

        roundTripEncoding(list, order1, decoder);

    }

    @Test
    public void testOrder1() throws IOException {
        IntArrayList list = new IntArrayList(list1);
        int numSymbols = 5;
        FastArithmeticCoderI order1 = new FastArithmeticCoderOrder1(numSymbols,list.size());
        FastArithmeticDecoderI decoder = new FastArithmeticDecoderOrder1(numSymbols,list.size());

        roundTripEncoding(list, order1, decoder);

    }

    @Test
       public void testPlus() throws IOException {
           IntArrayList list = new IntArrayList(list1);
           int numSymbols = 5;
           FastArithmeticCoderI order1 = new FastArithmeticCoderPlus(numSymbols,list.size());
           FastArithmeticDecoderI decoder = new FastArithmeticDecoderPlus(numSymbols, list.size());

           roundTripEncoding(list, order1, decoder);

       }

    private void roundTripEncoding(IntArrayList list, FastArithmeticCoderI order1, FastArithmeticDecoderI decoder) throws IOException {
        final FastByteArrayOutputStream arrayOutputStream = new FastByteArrayOutputStream();
        OutputBitStream out = new DebugOutputBitStream(new OutputBitStream(arrayOutputStream));
        for (int v : list) {
            order1.encode(v, out);
        }
        order1.flush(out);
        out.flush();

        final byte[] source = arrayOutputStream.array;
        byte[] buffer=new byte[source.length+10];
        System.arraycopy(source, 0, buffer, 0,source.length);
        System.out.printf("Encoded in %d bits%n", out.writtenBits());
        InputBitStream in = new DebugInputBitStream(new InputBitStream(buffer));

        for (int v : list) {

            final int x = decoder.decode(in);
            System.out.println("decoding: " + x);
            assertEquals(v, x);
        }
    }
}
