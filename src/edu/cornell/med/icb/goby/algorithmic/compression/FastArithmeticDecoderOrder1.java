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
import it.unimi.dsi.io.InputBitStream;

import java.io.IOException;

/**
 * An arithmetic decoder that models first order symbol dependencies in a naive way. See FastArithmeticDecoderOrder1
 *
 * @see FastArithmeticCoderOrder1
 *      Date: 1/15/12
 *      Time: 11:02 AM
 */
public final class FastArithmeticDecoderOrder1 implements FastArithmeticDecoderI {
    FastArithmeticDecoderI delegates[];
    private int numSymbols;
    private int previousSymbol;
    IntArrayList[] decodedLists;
    private boolean decoded = false;
    private int listSize;
    private int[] currentIndex;
    private int symbolRetrievalCount;
    private int[] lengths;

    public FastArithmeticDecoderOrder1(int numSymbols, int listSize) {
        delegates = new FastArithmeticDecoderI[numSymbols];
        decodedLists = new IntArrayList[numSymbols];
        lengths = new int[numSymbols];
        for (int i = 0; i < numSymbols; i++) {
            delegates[i] = new FastArithmeticDecoder(numSymbols);
            decodedLists[i] = new IntArrayList();
        }
        currentIndex = new int[numSymbols];
        this.numSymbols = numSymbols;
        this.listSize = listSize;
    }

    @Override
    public void reset() {
        previousSymbol = 0;
        decoded = false;
        int i = 0;
        for (FastArithmeticDecoderI delegate : delegates) {
            delegate.reset();
            decodedLists[i].clear();
            i += 1;
        }
    }

    @Override
    public int decode(InputBitStream ibs) throws IOException {
        if (!decoded) {
            for (int delegateIndex = 0; delegateIndex < numSymbols; delegateIndex++) {
                final int size = lengths[delegateIndex] = ibs.readNibble();
          //      System.out.printf("Reading %d symbols for predecessor=%d %n", size, delegateIndex);
                for (int i = 0; i < size; i++) {
                    int x = delegates[delegateIndex].decode(ibs);
                    decodedLists[delegateIndex].add(x);
                }
                //        System.out.printf("%n1. delegate %d is now positioned at %d %n", delegateIndex, ibs.readBits());
                reposition(ibs, delegateIndex);

              //  System.out.printf("2. delegate %d is now positioned at %d %n", delegateIndex, ibs.readBits());
            }
            decoded = true;
        }
        assert symbolRetrievalCount < listSize : "You cannot retrieve more than listSize symbols";
        final int symbol = decodedLists[previousSymbol].get(currentIndex[previousSymbol]++);
        previousSymbol = symbol;
        symbolRetrievalCount += 1;
        return symbol;
    }

    @Override
    public void flush(InputBitStream ibs) throws IOException {
        throw new UnsupportedOperationException("flush is not supported by this implementation.");

        /*    for (FastArithmeticDecoderI delegate : delegates) {
         delegate.flush(ibs);
     }
     ibs.flush();   */
    }


    public void flush(InputBitStream input, int delegateIndex) throws IOException {
        delegates[delegateIndex].flush(input);
    }

    @Override
    public long getWindow() {
        throw new UnsupportedOperationException("getWindow is not supported by this implementation.");
    }

    @Override
    public void reposition(InputBitStream input) throws IOException {
        throw new UnsupportedOperationException("reposition is not supported by this implementation.");

        /*for (FastArithmeticDecoderI delegate : delegates) {
            delegate.reposition(input);
        } */
    }

    public void reposition(InputBitStream input, int delegateIndex) throws IOException {
        flush(input, delegateIndex);
        final long readBits = input.readBits();

        final long position = readBits - FastArithmeticDecoder.BITS;
        //       System.out.printf("readBits= %d%n",  input.readBits());
        input.flush();
        if (position >= 0) {
            input.position(position);
            input.readBits(position);
        }
    }

}
