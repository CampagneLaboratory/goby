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

/*
 * MG4J: Managing Gigabytes for Java (big)
 *
 * Copyright (C) 2002-2011 Sebastiano Vigna
 *
 *  This library is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This library is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
 *  for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 */

import it.unimi.dsi.io.OutputBitStream;

import java.io.IOException;

/**
 * An arithmetic coder adapted from MG4J's implementation. This implementation adds
 * methods to reset the coder to support encoding multiple messages interleaved with
 * other streams. This makes it possible to continually refine the symbol frequencies
 * over multiple messages (which would not happen if we created a new instance of the
 * coder for each message, and keep synchronization in the output stream. The class is
 * called FastArithmeticCoder for symmetry with the FastArithmeticDecoder, but is not
 * faster than the MG4J class.
 *
 * @author Sebastiano Vigna
 * @author Fabien Campagne
 * @see it.unimi.dsi.mg4j.io.ArithmeticDecoder (MG4J).
 * @since 0.1
 */

final public class FastArithmeticCoder implements FastArithmeticCoderI {
    /**
     * Number of bits used by the coder.
     */
    public final static int BITS = 63;

    /**
     * Bit-level representation of 1/2.
     */
    private final static long HALF = 1L << (BITS - 1);

    /**
     * Bit-level representation of 1/4.
     */
    private final static long QUARTER = 1L << (BITS - 2);

    /**
     * Cumulative counts for all symbols.
     */
    private int[] cumCount;

    /**
     * Total count.
     */
    private int total;

    /**
     * Number of symbols.
     */
    private int n;

    /**
     * Current base of the range.
     */
    private long low = 0;

    /**
     * Current width of the range.
     */
    private long range = HALF;

    /**
     * Bits held from output.
     */
    private long outstandingBits = 0;

    /**
     * The first bit is always 0, so we do not output it.
     */
    private boolean firstBit = true;

    @Override
    public void reset() {
        range = HALF;
        outstandingBits = 0;
        low = 0;
        firstBit = true;
        /* cumCount = new int[ n + 1 ];
          for ( int i = 0; i < n; i++ )
              incrementCount( i ); // Initially, everything is equiprobable.
          total = n;*/
    }

    /**
     * Creates a new coder.
     *
     * @param n number of symbols used by the coder.
     */

    public FastArithmeticCoder(final int n) {
        if (n < 1)
            throw new IllegalArgumentException("You cannot use " + n + " symbols.");
        this.n = n;
        cumCount = new int[n + 1];
        for (int i = 0; i < n; i++)
            incrementCount(i); // Initially, everything is equiprobable.
        total = n;
    }


    /* The following methods implement a Fenwick tree. */

    private void incrementCount(final int xArg) {
        int x=xArg;
        x++;

        while (x <= n) {
            cumCount[x]++;
            x += x & -x; // By chance, this gives the right next index 8^).
        }
    }

    private int getCount(final int xArg) {

        int x = xArg;
        int c = 0;
        while (x != 0) {
            c += cumCount[x];
            x = x & x - 1; // This cancels out the least nonzero bit.
        }

        return c;
    }


    /**
     * Writes a bit and all outstanding bits.
     * <p/>
     * <P>First the given bit is output. Then, {@link #outstandingBits} opposite bits are output.
     * The first overall bit (which is always 0) is never output.
     *
     * @param bit a bit.
     * @param obs the output stream.
     * @return the number of bits written.
     */

    private int emit(final int bit, final OutputBitStream obs) throws IOException {

        if (firstBit) {
            firstBit = false;
            return 0;
        }

        int l = obs.writeBit(bit);
        while (outstandingBits-- != 0)
            l += obs.writeBit(1 - bit);
        outstandingBits = 0;
        return l;
    }

    /**
     * Encodes a symbol.
     *
     * @param x   a bit.
     * @param obs the output stream.
     * @return the number of bits written (note that it can be 0, as arithmetic compression can
     *         encode a symbol in a fraction of a bit).
     * @throws IOException if <code>obs</code> does.
     */

    @Override
    public int encode(final int x, final OutputBitStream obs) throws IOException {
        if (x < 0)
            throw new IllegalArgumentException("You cannot encode a negative symbol.");
        if (x >= n)
            throw new IllegalArgumentException("You cannot encode " + x + ": you have only " + n + " symbols.");

        final long r = range / total;
        final int lowCount = getCount(x);
        final int highCount = getCount(x + 1);

        low += r * lowCount;

        if (x != n - 1)
            range = r * (highCount - lowCount);
        else
            range -= r * lowCount;

        incrementCount(x);
        total++;

        int l = 0;

        while (range <= QUARTER) {
            if (low >= HALF) {
                l += emit(1, obs);
                low -= HALF;
            } else if (range + low <= HALF) {
                l += emit(0, obs);
            } else {
                low -= QUARTER;
                outstandingBits++;
            }
            range <<= 1;
            low <<= 1;
        }

        return l;
    }

    /**
     * Flushes the last bits.
     * <p/>
     * <P>This method must be called when coding is over. It guarantees that enough
     * bits are output to make the decoding of the last symbol nonambiguous, whichever
     * bits follow in the stream.
     *
     * @param obs the output stream.
     * @return the number of bits written.
     * @throws IOException if <code>obs</code> does.
     */

    @Override
    public int flush(final OutputBitStream obs) throws IOException {
        int nbits, i, l = 0;
        long roundup, bits = 0, value;

        for (nbits = 1; nbits <= BITS; nbits++) {
            roundup = (1L << (BITS - nbits)) - 1;
            bits = (low + roundup) >>> (BITS - nbits);
            value = bits << (BITS - nbits);

            if (low <= value && (value + roundup <= low + (range - 1) || value + roundup >= 0 && low + (range - 1) < 0) // This handles overflows onto the most significant bit.
                    )
                break;
        }

        for (i = 1; i <= nbits; i++)
            l += emit((int) ((bits >>> (nbits - i)) & 1), obs);

        return l;
    }

}
