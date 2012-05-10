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
import org.apache.tools.bzip2.CBZip2InputStream;
import org.apache.tools.bzip2.CBZip2OutputStream;

import java.io.*;

/**
 * A BZIP2 Chunk coder. Simply bzip2 the protocol buffer collection.
 *
 * @author Fabien Campagne
 *         Date: 3/8/12
 *         Time: 8:31 AM
 */
public class BZip2ChunkCodec implements ChunkCodec {

    private ProtobuffCollectionHandler parser;

    @Override
    public void setHandler(final ProtobuffCollectionHandler parser) {
        this.parser = parser;
    }

    private final byte[] bytes=new byte[7];

    @Override
    public boolean validate(DataInputStream input) {
        try {
            final int length = 4 + 3;    // size 4 bytes + magic number 1F 8B 08

            if (input.read(bytes, 0, length) != length) {
                return false;
            } else {
                return bytes[4] == (byte)0x42 && bytes[5] == (byte)0x5A && bytes[6] == (byte)0x68;
            }
        } catch (IOException e) {
            return false;
        }
    }

    @Override
    public String name() {
        return "bzip2";
    }

    @Override
    public byte registrationCode() {
        return REGISTRATION_CODE;
    }

    public static final byte REGISTRATION_CODE = -3;

    @Override
    public ByteArrayOutputStream encode(final Message readCollection) throws IOException {
        final ByteArrayOutputStream byteBuffer = new ByteArrayOutputStream(10000);

        final OutputStream bZip2OutputStream = new CBZip2OutputStream(byteBuffer);
        readCollection.writeTo(bZip2OutputStream);
        bZip2OutputStream.flush();
        bZip2OutputStream.close();
        return byteBuffer;

    }

    @Override
    public Message decode(final byte[] bytes) throws IOException {
        final CBZip2InputStream uncompressStream = new CBZip2InputStream(new ByteArrayInputStream(bytes));
        try {
            return parser.parse(uncompressStream);
        } finally {
            uncompressStream.close();
        }


    }

    @Override
    public int getSuggestedChunkSize() {
        return 20000;
    }

}
