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

package edu.cornell.med.icb.goby.compression;

import com.google.protobuf.Message;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

/**
 *  A null chunk codec that does not write anything. Only useful to measure time taken by other processes for benchmarks.
 *  Please note that you will be unable to recover anything compressed with this codec, it all goes to /dev/null!
 *  @author Fabien Campagne
 *         Date: 3/18/12
 *         Time: 10:13 AM
 */
public class NullChunkCodec implements ChunkCodec {

    private ProtobuffCollectionHandler parser;

    @Override
    public void setHandler(final ProtobuffCollectionHandler parser) {
        this.parser = parser;
    }

    @Override
    public String name() {
        return "null";
    }

    @Override
    public byte registrationCode() {
        return  REGISTRATION_CODE;
    }

    public static final byte REGISTRATION_CODE = (byte) -4;

    @Override
    public ByteArrayOutputStream encode(final Message readCollection) throws IOException {
        return new ByteArrayOutputStream(0);

    }

    @Override
    public Message decode(final byte[] bytes) throws IOException {
        throw new UnsupportedOperationException("The null codec cannot retrieve any data (it all went to /dev/null!), so decode is not supported.");

    }

}
