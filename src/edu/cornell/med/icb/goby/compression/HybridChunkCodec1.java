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
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.*;

/**
 * A codec that writes highly compressed data in one pool, but keeps protobuf message in a separate pool.
 * @author Fabien Campagne
 *         Date: 3/3/12
 *         Time: 2:35 PM
 */
public class HybridChunkCodec1 implements ChunkCodec {
    private boolean debug = false;

    @Override
    public String name() {
        return "hybrid-1";
    }

    public static final byte REGISTRATION_CODE = -2;
    private ProtobuffCollectionHandler handler;
    private final GZipChunkCodec gzipCodec = new GZipChunkCodec();

    @Override
    public byte registrationCode() {
        return REGISTRATION_CODE;
    }

    /**
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(HybridChunkCodec1.class);
    private int chunkIndex = 0;

    @Override
    public ByteArrayOutputStream encode(final Message readCollection) throws IOException {
        if (readCollection == null) {
            return null;
        }
        final ByteArrayOutputStream result = new ByteArrayOutputStream();
        final DataOutputStream completeChunkData = new DataOutputStream(result);

        final ByteArrayOutputStream compressedBits = new ByteArrayOutputStream();
        final Message reducedProtoBuff = handler.compressCollection(readCollection, compressedBits);

        final int compressedBitSize = compressedBits.size();
        completeChunkData.writeInt(compressedBitSize);
        completeChunkData.write(compressedBits.toByteArray());

        final ByteArrayOutputStream out = gzipCodec.encode(reducedProtoBuff);

        final byte[] gzipBytes = out.toByteArray();
        final int gzipBytesSize = gzipBytes.length;
        completeChunkData.write(gzipBytes);
        completeChunkData.flush();
        if (debug && chunkIndex % 100 == 0) {

            //TODO remove compression of original collection. Only useful for stat collection
            int originalGzipSize = gzipCodec.encode(readCollection).toByteArray().length;

            final int gain = originalGzipSize - (gzipBytesSize + compressedBitSize);
            LOG.info(String.format("compressed size=%d gzip size=%d (original gzip=%d) percent compressed/(compressed+gzip) %g gain=%d, %g%% ",
                    compressedBitSize, gzipBytesSize, originalGzipSize,
                    100d * ((double) compressedBitSize) / (compressedBitSize + gzipBytesSize),
                    gain, gain * 100d / originalGzipSize));

        }
        chunkIndex++;
        return result;
    }

    @Override
    public Message decode(final byte[] bytes) throws IOException {
        final DataInputStream completeChunkData = new DataInputStream(new ByteArrayInputStream(bytes));
        final int compressedSize = completeChunkData.readInt();
        final byte[] compressedBytes = new byte[compressedSize];
        final int read = completeChunkData.read(compressedBytes, 0, compressedSize);
        assert read == compressedSize : "read size must match recorded size.";
        final int bytesLeft = bytes.length - 4 - compressedSize;
        final byte[] leftOver = new byte[bytesLeft];
        // 4 is the number of bytes to encode the length of the compressed chunk.
        System.arraycopy(bytes, 4 + compressedSize, leftOver, 0, bytesLeft);
        final Message reducedProtoBuff = gzipCodec.decode(leftOver);
        return handler.decompressCollection(reducedProtoBuff, compressedBytes);
    }

    @Override
    public void setHandler(final ProtobuffCollectionHandler handler) {
        this.handler = handler;
        gzipCodec.setHandler(handler);
    }


}
