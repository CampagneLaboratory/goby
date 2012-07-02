/*
 *
 * Copyright (C) 2002-2011 Sebastiano Vigna, 2012 Fabien Campagne
 *
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

package edu.cornell.med.icb.goby.algorithmic.compression;

import it.unimi.dsi.bits.Fast;
import it.unimi.dsi.io.InputBitStream;

import java.io.IOException;

/**
 * An arithmetic decoder that works well with large alphabets of symbols (>100 symbols). This class is copied and
 * modified from the MG4J distribution. It was changed to use binary search in the decode method where MG4J used a
 * linear search. The modification provides large increase in performance when the alphabet is large (>100).
 *
 * @author Sebastiano Vigna (original class from MG4J)
 * @author Fabien Campagne (modification to improve scalability for large alphabets)
 * @see it.unimi.dsi.mg4j.io.ArithmeticDecoder (MG4J)
 *      Date: 1/15/12
 *      Time: 11:02 AM
 */
public final class FastArithmeticDecoder implements FastArithmeticDecoderI {
    /**
     * Number of bits used by the decoder.
     */
    public static final int BITS = 63;

    /**
     * Bit-level representation of 1/2.
     */
    private static final long HALF = 1L << (BITS - 1);

    /**
     * Bit-level representation of 1/4.
     */
    private static final long QUARTER = 1L << (BITS - 2);

    /**
     * Cumulative counts for all symbols.
     */
    private final int count[];

    /**
     * Total count.
     */
    private int total;

    /**
     * Number of symbols.
     */
    private final int n;

    /**
     * Current width of the range.
     */
    private long range = HALF;

    /**
     * Current bits being decoded.
     */
    private long buffer = -1;

    /**
     * Current window on the bit stream.
     */
    private long window = 0;
    private boolean useBinarySearch = true;
    private int maxIndexNonZeroFrequency = Integer.MIN_VALUE;
    private int lastCount;
    private int lastIndex=-1;

    /**
     * Resets the decoder before decoding a new message. The method prepares the coder for the first call
     * to encode. Frequencies of symbols are unchanged after calling this method.
     */
    @Override
    public void reset() {
        window = 0;
        buffer = -1;
        range = HALF;

    }

    /**
     * Creates a new decoder.
     *
     * @param n number of symbols used by the decoder.
     */

    public FastArithmeticDecoder(final int n) {
        if (n < 1)
            throw new IllegalArgumentException("You cannot use " + n + " symbols.");
        this.n = n;
        count = new int[n + 1];
        for (int i = 0; i < n; i++) {
            incrementCount(i); // Initially, everything is equiprobable.
        }
        total = n;
    }


    /* The following methods implement a Fenwick tree. */

    private void incrementCount(int x) {

        x++;
        while (x <= n) {
            count[x]++;
            x += x & -x; // By chance, this gives the right next index 8^).
        }

    }


    private int getCount(int x) {
        if (x == lastIndex) {
            // used the last cached count for x, as calculated previously in findXBinary
            return lastCount;
        }
        int c = 0;
        final int xArg = x;
        while (x != 0) {
            c += count[x];
            x = x & x - 1; // This cancels out the least nonzero bit.
        }

        return c;
    }

    private int findXBinary(final int x, final int start, final int end) {
        if (end == start) {
       //     System.out.printf("Returning from findX with x=%d %n", start - 1);
            return start - 1;
        }
        final int middle = (start + end) / 2;
        final int count = getCount(middle);
       // System.out.printf("start=%d middle=%d end=%d count=%d %n", start, middle, end, count);
        if (x < count) {

            return findXBinary(x, start, middle);
        }
        if (x > count) {
            return findXBinary(x, middle + 1, end);
        }
     //   System.out.printf("Returning from findX with x=%d %n", middle);
        lastCount = count;
        lastIndex = middle;
        return middle;
    }


    /**
     * Decodes a symbol.
     *
     * @param ibs the input stream.
     * @return the next symbol encoded.
     * @throws IOException if <code>ibs</code> does.
     */

    @Override
    public int decode(final InputBitStream ibs) throws IOException {

        if (buffer == -1) {
            window = buffer = ibs.readLong(BITS - 1); // The first output bit is always 0 and is not output.
        }
        final long r = range / total;
        int x = (int) (buffer / r);
        if (total - 1 < x) {
            x = total - 1;
        }
        final int v = x;

        // We replace the above code with a binary search algorithm (O(log N) for alphabets of size N).

        x = findXBinary(v, 1, n);

        final int lowCount = getCount(x);
        final int highCount = getCount(x + 1);
        // invalidate the last count cache:
        lastIndex=-1;
        final long l = r * lowCount;
        buffer -= l;

        if (x != n - 1) {
            range = r * (highCount - lowCount);
        } else {
            range -= l;
        }
        incrementCount(x);
        total++;

        while (range <= QUARTER) {
            buffer <<= 1;
            range <<= 1;
            window <<= 1;
            if (ibs.readBit() != 0) {
                buffer += 1;
                window += 1;
            }
        }

        return x;
    }


    /**
     * Flushes (reads) the disambiguating bits.
     * <p/>
     * <P>This method must be called when all symbols have been decoded.  After
     * the call, exactly {@link #BITS} excess bits will be present in the
     * current window of the decoder. You can get them using {@link #getWindow()};
     * usually you will then unget them to the bit stream.
     *
     * @param ibs the input stream.
     * @throws IOException if <code>ibs</code> does.
     */

    @Override
    public void flush(final InputBitStream ibs) throws IOException {
        int nbits, i;
        long roundup, bits, value, low;

        low = ((window & (HALF - 1)) + HALF) - buffer;

        for (nbits = 1; nbits <= BITS; nbits++) {
            roundup = (1L << (BITS - nbits)) - 1;
            bits = (low + roundup) >>> (BITS - nbits);
            value = bits << (BITS - nbits);

            if (low <= value && (value + roundup <= low + (range - 1) || value + roundup >= 0 && low + (range - 1) < 0) // This handles overflows onto the most significant bit.
                    )
                break;
        }

        for (i = 1; i <= nbits; i++) {
            window <<= 1;
            window |= ibs.readBit();
        }

        // System.out.println("flushed nbits:"+nbits);


    }

    /**
     * Returns the current bit stream window.
     *
     * @return the current bit stream window in the lower {@link #BITS} bits.
     */

    @Override
    public long getWindow() {
        return window & ((HALF << 1) - 1);
    }

    /**
     * Reposition the input stream. Repositioning is needed after decoding to make it possible to read from the
     * stream with with decoder. Calling this method flushes the decoder and repositions the input stream
     * at the exact end of the flush bits.
     *
     * @param input
     * @throws IOException
     */
    @Override
    public void reposition(InputBitStream input) throws IOException {
        flush(input);
        final long position = input.readBits() - FastArithmeticDecoder.BITS;
        //       System.out.printf("readBits= %d%n",  input.readBits());
        input.flush();
        input.position(position);
        input.readBits(position);
    }
}
