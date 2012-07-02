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
 * An arithmetic coder that models first order symbol dependencies in a naive way. The probability of symbol
 * x is modeled as a function of the symbol encoded immediately before x. All pairs of possible symbols are
 * considered in this naive implementation. Note the drawback that if some symbol pairs just never occur in
 * a given input, they are represented in the n^2 fenwick tree and will degrade the probability estimates
 * for the pairs that do occur.
 *
 * @author Fabien Campagne
 * @see FastArithmeticDecoderOrder1
 * @since 0.1
 */

final public class FastArithmeticCoderOrder1 implements FastArithmeticCoderI {
    FastArithmeticCoderI delegates[];
    private int numSymbols;
    private int previousSymbol;

    private OutputBitStream[] outs;
       private FastByteArrayOutputStream[] arrays;
    /** Number of symbols written that follow previousSymbol:
     *
     */
    private int[] lengths;

    public FastArithmeticCoderOrder1(final int numSymbols) {
        delegates = new FastArithmeticCoderI[numSymbols];
        outs=new OutputBitStream[numSymbols];
        arrays=new FastByteArrayOutputStream[numSymbols];
        lengths=new int[numSymbols];
        for (int i = 0; i < numSymbols; i++) {
            delegates[i] = new FastArithmeticCoder(numSymbols);
            arrays[i]=new FastByteArrayOutputStream();
            outs[i]=new OutputBitStream(arrays[i]);
        }
        this.numSymbols = numSymbols;

    }



    @Override
    public void reset() {
        previousSymbol=0;
        for (final FastArithmeticCoderI delegate : delegates) {
            delegate.reset();
        }
    }

    @Override
    public int encode(final int x, final OutputBitStream obs) throws IOException {
        lengths[previousSymbol]++;
        final int result = delegates[previousSymbol].encode(x, outs[previousSymbol]);
        previousSymbol = x;
        return result;
    }

    @Override
    public int flush(final OutputBitStream obs) throws IOException {
        int bitsWritten = 0;

        int i=0;
        for (final FastArithmeticCoderI delegate : delegates) {
            int delegateBitsWritten = 0;
            int bitsFlushed=delegate.flush(outs[i]);

           // System.out.printf("bits flushed: %d %n",bitsFlushed);
            final long length = outs[i].writtenBits();
            outs[i].flush();
            // write the number of symbols output for a given predecessorSymbol:
            delegateBitsWritten+=obs.writeNibble(lengths[i]);
            delegateBitsWritten+=obs.write(arrays[i].array, length);
         //   System.out.printf("delegate %d wrote %d bits %n",i, delegateBitsWritten);
            i++;
            bitsWritten+=delegateBitsWritten;
        }
    //    System.out.printf("all delegates wrote %d bits %n",bitsWritten);
       // obs.flush();
        return bitsWritten;
    }
}
