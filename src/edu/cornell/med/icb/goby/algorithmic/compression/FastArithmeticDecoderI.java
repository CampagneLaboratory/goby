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

import it.unimi.dsi.io.InputBitStream;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: 6/30/12
 *         Time: 2:53 PM
 */
public interface FastArithmeticDecoderI {
    /**
     * Resets the decoder before decoding a new message. The method prepares the coder for the first call
     * to encode. Frequencies of symbols are unchanged after calling this method.
     */
    void reset();

    /**
     * Decodes a symbol.
     *
     * @param ibs the input stream.
     * @return the next symbol encoded.
     * @throws java.io.IOException if <code>ibs</code> does.
     */

    int decode(InputBitStream ibs) throws IOException;

    /**
     * Flushes (reads) the disambiguating bits.
     * <p/>
     * <P>This method must be called when all symbols have been decoded.  After
     * the call, exactly {@link #BITS} excess bits will be present in the
     * current window of the decoder. You can get them using {@link #getWindow()};
     * usually you will then unget them to the bit stream.
     *
     * @param ibs the input stream.
     *            @return  the number of bits flushed
     * @throws java.io.IOException if <code>ibs</code> does.
     */

    void flush(InputBitStream ibs) throws IOException;

    /**
     * Returns the current bit stream window.
     *
     * @return the current bit stream window in the lower {@link #BITS} bits.
     */

    long getWindow();

    /**
     * Reposition the input stream. Repositioning is needed after decoding to make it possible to read from the
     * stream with with another encoder. Calling this method flushes the decoder and repositions the input stream
     * at the exact end of the flush bits.
     *
     * @param input
     * @throws java.io.IOException
     */
    void reposition(InputBitStream input) throws IOException;
}
