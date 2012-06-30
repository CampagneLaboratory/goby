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

import it.unimi.dsi.fastutil.io.FastByteArrayOutputStream;
import it.unimi.dsi.io.OutputBitStream;

import java.io.IOException;

/**
 * An arithmetic coder that models symbol probabilities differently when the previous symbol was the
 * most abundant seen so far in the input.
 *
 * @author Fabien Campagne
 * @see edu.cornell.med.icb.goby.algorithmic.compression.FastArithmeticDecoderPlus
 * @since 0.1
 */

final public class FastArithmeticCoderPlus implements FastArithmeticCoderI {
    FastArithmeticCoder delegates[];
    private int numSymbols;
    private int previousSymbol;

    /**
     * The exact number of symbols that will be encoded with this coder, before the next call to reset().
     */
    private int listSize;
    private OutputBitStream[] outs;
    private FastByteArrayOutputStream[] arrays;
    /**
     * Number of symbols written that follow previousSymbol:
     */
    private int[] lengths;
    private int mostAbundantCount;
    private int mostAbundantSymbol = -1;
    private static final int MOST_ABUNDANT_INDEX = 0;
    private int[] counts;

    public FastArithmeticCoderPlus(final int numSymbols, final int listSize) {
        final int numCoders = 2;
        delegates = new FastArithmeticCoder[numCoders];
        outs = new OutputBitStream[numCoders];
        arrays = new FastByteArrayOutputStream[numCoders];
        lengths = new int[numCoders];
        counts = new int[numSymbols];
        for (int i = 0; i < numCoders; i++) {
            delegates[i] = new FastArithmeticCoder(numSymbols);
            arrays[i] = new FastByteArrayOutputStream();
            outs[i] = new OutputBitStream(arrays[i]);
        }
        this.numSymbols = numSymbols;
        this.listSize = listSize;
    }


    @Override
    public void reset() {
        previousSymbol = 0;
        for (final FastArithmeticCoderI delegate : delegates) {
            delegate.reset();
        }
    }

    @Override
    public int encode(final int x, final OutputBitStream obs) throws IOException {
        int zeroOrderIndex = 1;
        final int result;
        final int previousCount = counts[previousSymbol];
        if (previousCount > mostAbundantCount || mostAbundantSymbol == previousSymbol) {
            mostAbundantSymbol = previousSymbol;
            mostAbundantCount = previousCount;
            // x is the single most abundant symbol encountered so far in the input. Switch to order 1 modeling:
            lengths[MOST_ABUNDANT_INDEX]++;
            result = delegates[MOST_ABUNDANT_INDEX].encode(x, outs[MOST_ABUNDANT_INDEX]);

        } else {
            result = delegates[zeroOrderIndex].encode(x, outs[zeroOrderIndex]);
            lengths[zeroOrderIndex]++;
        }
        previousSymbol = x;
        counts[x]++;

        return result;
    }

    @Override
    public int flush(final OutputBitStream obs) throws IOException {
        int bitsWritten = 0;

        int i = 0;
        for (final FastArithmeticCoderI delegate : delegates) {
            int delegateBitsWritten = 0;
            int bitsFlushed = delegate.flush(outs[i]);

           // System.out.printf("bits flushed: %d %n", bitsFlushed);
            final long length = outs[i].writtenBits();
            outs[i].flush();
            // write the number of symbols output for a given predecessorSymbol:
            delegateBitsWritten += obs.writeNibble(lengths[i]);
            delegateBitsWritten += obs.write(arrays[i].array, length);
      //      System.out.printf("delegate %d wrote %d bits %n", i, delegateBitsWritten);
            i++;
            bitsWritten += delegateBitsWritten;
        }
        //    System.out.printf("all delegates wrote %d bits %n",bitsWritten);
        // obs.flush();
        return bitsWritten;
    }
}
