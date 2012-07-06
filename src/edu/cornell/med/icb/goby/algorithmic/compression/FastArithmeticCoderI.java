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

import it.unimi.dsi.io.OutputBitStream;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: 6/30/12
 *         Time: 2:52 PM
 */
public interface FastArithmeticCoderI {
    void reset();

    /** Encodes a symbol.
     *
     * @param x a bit.
     * @param obs the output stream.
     * @return the number of bits written (note that it can be 0, as arithmetic compression can
     * encode a symbol in a fraction of a bit).
     * @throws java.io.IOException if <code>obs</code> does.
     */

    int encode(int x, OutputBitStream obs) throws IOException;

    /** Flushes the last bits.
     *
     * <P>This method must be called when coding is over. It guarantees that enough
     * bits are output to make the decoding of the last symbol nonambiguous, whichever
     * bits follow in the stream.
     *
     * @param obs the output stream.
     * @return the number of bits written.
     * @throws java.io.IOException if <code>obs</code> does.
     */

    int flush(OutputBitStream obs) throws IOException;
}
