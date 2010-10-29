/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.reads;

import java.util.Arrays;
import java.util.zip.CRC32;
import java.util.zip.Checksum;


/**
 * Compressed representation of the sequence of a read.
 *
 * @author Fabien Campagne
 *         Date: Jun 3, 2009
 *         Time: 2:12:28 PM
 */
public class CompressedRead {
    private final byte[] data;

    private static final Checksum DIGESTER = new CRC32();
    private final int hashCode;

    /**
     * Read index is not part of the hashCode and not considerd by the equals method. It is
     * only here to help map back to the position of the read in the inoput file.
     */
    public int readIndex;

    @Override
    public final int hashCode() {
        return hashCode;
    }

    @Override
    public final boolean equals(final Object o) {
        if (!(o instanceof CompressedRead)) {
            return false;
        }
        final CompressedRead other = (CompressedRead) o;

        return Arrays.equals(data, other.data);

    }

    public CompressedRead(final byte[] byteBuffer) {
        this.data = byteBuffer;
        hashCode = hashCode(this.data);
    }

    public static int hashCode(final byte[] data) {
        synchronized (DIGESTER) {
            final long digest = hashCodeLong(data);
            return ((32 >> digest)) ^ ((int) (digest));
        }
    }

    public static long hashCodeLong(final byte[] data) {
        synchronized (DIGESTER) {
            DIGESTER.reset();
            DIGESTER.update(data, 0, data.length);
            return DIGESTER.getValue();
        }
    }

    public static long hashCodeLongNotSynchronized(final byte[] data) {
        DIGESTER.reset();
        DIGESTER.update(data, 0, data.length);
        return DIGESTER.getValue();
    }

    public final long hashCodeLong() {
        return hashCodeLong(this.data);
    }
}
