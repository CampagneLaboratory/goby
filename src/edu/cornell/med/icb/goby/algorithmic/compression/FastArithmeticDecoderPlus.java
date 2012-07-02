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
 * An arithmetic decoder that decodes streams generated with FastArithmeticCoderPlus.
 *
 * @see edu.cornell.med.icb.goby.algorithmic.compression.FastArithmeticCoderPlus
 *      Date: 1/15/12
 *      Time: 11:02 AM
 */
public final class FastArithmeticDecoderPlus implements FastArithmeticDecoderI {
    private final FastArithmeticDecoder delegates[];
    private int previousSymbol;
    private final IntArrayList[] decodedLists;
    private boolean decoded = false;
    private final int[] currentIndex;
    private int symbolRetrievalCount;
    private final int[] lengths;
    private final int numCoders;
    private int mostAbundantCount;
    private int mostAbundantSymbol = -1;
    private final int[] counts;

    public FastArithmeticDecoderPlus(final int numSymbols) {
        numCoders = 2;
        delegates = new FastArithmeticDecoder[numCoders];
        decodedLists = new IntArrayList[numCoders];
        lengths = new int[numCoders];
        for (int i = 0; i < numCoders; i++) {
            delegates[i] = new FastArithmeticDecoder(numSymbols);
            decodedLists[i] = new IntArrayList();
        }
        counts = new int[numSymbols];
        currentIndex = new int[numCoders];

    }

    @Override
    public void reset() {
        previousSymbol = 0;
        decoded = false;
        int i = 0;
        for (final FastArithmeticDecoderI delegate : delegates) {
            delegate.reset();
            decodedLists[i].clear();
            i += 1;
        }
    }

    @Override
    public int decode(final InputBitStream ibs) throws IOException {
        if (!decoded) {
            for (int delegateIndex = 0; delegateIndex < numCoders; delegateIndex++) {
                final int size = ibs.readNibble();
                final FastArithmeticDecoder delegate = delegates[delegateIndex];
                final IntArrayList decodedList = decodedLists[delegateIndex];
                if (size > 0) {
                    //      System.out.printf("Reading %d symbols for predecessor=%d %n", size, delegateIndex);
                    for (int i = 0; i < size; i++) {
                        final int x = delegate.decode(ibs);
                        decodedList.add(x);
                    }
                    //        System.out.printf("%n1. delegate %d is now positioned at %d %n", delegateIndex, ibs.readBits());
                    reposition(ibs, delegateIndex);
                }
                //  System.out.printf("2. delegate %d is now positioned at %d %n", delegateIndex, ibs.readBits());
            }
            decoded = true;
        }
        //  assert symbolRetrievalCount <= listSize : "You cannot retrieve more than listSize symbols";
        final int delegateIndex = previousSymbol == mostAbundantSymbol ? 0 : 1;
        final int symbol = decodedLists[delegateIndex].getInt(currentIndex[delegateIndex]++);
        previousSymbol = symbol;
        counts[previousSymbol] += 1;

        final int previousCount = counts[previousSymbol];
        if (mostAbundantSymbol == previousSymbol || previousCount > mostAbundantCount) {
            mostAbundantSymbol = previousSymbol;
            mostAbundantCount = previousCount;
        }
        return symbol;
    }

    @Override
    public void flush(final InputBitStream ibs) throws IOException {
        throw new UnsupportedOperationException("flush is not supported by this implementation.");

        /*    for (FastArithmeticDecoderI delegate : delegates) {
         delegate.flush(ibs);
     }
     ibs.flush();   */
    }


    public void flush(final InputBitStream input, final int delegateIndex) throws IOException {
        delegates[delegateIndex].flush(input);
    }

    @Override
    public long getWindow() {
        throw new UnsupportedOperationException("getWindow is not supported by this implementation.");
    }

    @Override
    public void reposition(final InputBitStream input) throws IOException {
        // we have already repositioned.
    }

    public void reposition(final InputBitStream input, final int delegateIndex) throws IOException {
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
