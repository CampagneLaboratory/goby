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
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

/**
 *  A Chunk coder that does not compress. Useful only as a baseline to determine processing times
 *  without compression.
 *  @author Fabien Campagne
 *         Date: 3/9/12
 *         Time: 9:10 AM
 */
public class NotCompressedChunkCodec implements ChunkCodec {

    private ProtobuffCollectionHandler parser;

    @Override
    public void setHandler(final ProtobuffCollectionHandler parser) {
        this.parser = parser;
    }

    @Override
    public boolean validate(FastBufferedInputStream input) {
        return true;
    }

    @Override
    public String name() {
        return "no-compression";
    }

    @Override
    public byte registrationCode() {
        return  REGISTRATION_CODE;
    }

    public static final byte REGISTRATION_CODE = -4;

    @Override
    public ByteArrayOutputStream encode(final Message readCollection) throws IOException {
        final ByteArrayOutputStream byteBuffer = new ByteArrayOutputStream(10000);
        readCollection.writeTo(byteBuffer);
        byteBuffer.flush();
        byteBuffer.close();
        return byteBuffer;

    }

    @Override
    public Message decode(final byte[] bytes) throws IOException {
        final ByteArrayInputStream uncompressStream = new ByteArrayInputStream(bytes);
        try {
            return parser.parse(uncompressStream);
        } finally {
            uncompressStream.close();
        }


    }

    @Override
    public int getSuggestedChunkSize() {
        return 10000;
    }

}
