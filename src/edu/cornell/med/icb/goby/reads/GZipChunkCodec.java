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

package edu.cornell.med.icb.goby.reads;

import com.google.protobuf.GeneratedMessage;
import com.google.protobuf.Message;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 *  The original Goby Gzip Chunk coder. Simply GZips the protocol buffer collection.
 *  @author Fabien Campagne
 *         Date: 3/3/12
 *         Time: 10:30 AM
 */
public class GZipChunkCodec implements ChunkCodec {

    private ProtobuffCollectionParser parser;

    @Override
    public void setParser(final ProtobuffCollectionParser parser) {
        this.parser = parser;
    }

    @Override
    public String name() {
        return "gzip";
    }

    @Override
    public byte registrationCode() {
        return (byte) 0xFF;
    }

    public static final byte REGISTRATION_CODE = (byte) 0xFF;

    @Override
    public ByteArrayOutputStream encode(final Message readCollection) throws IOException {
        final ByteArrayOutputStream byteBuffer = new ByteArrayOutputStream(10000);

        final OutputStream gzipOutputStream = new GZIPOutputStream(byteBuffer);
        readCollection.writeTo(gzipOutputStream);
        gzipOutputStream.flush();
        gzipOutputStream.close();
        return byteBuffer;

    }

    @Override
    public GeneratedMessage decode(final byte[] bytes) throws IOException {
        final GZIPInputStream uncompressStream = new GZIPInputStream(new ByteArrayInputStream(bytes));
        try {
            return parser.parse(uncompressStream);
        } finally {
            uncompressStream.close();
        }


    }

}
