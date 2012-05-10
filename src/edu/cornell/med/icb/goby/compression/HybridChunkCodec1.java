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
import java.util.zip.CRC32;

/**
 * A codec that writes highly compressed data in one pool, and keeps left-over uncompressed protobuf messages in a separate pool.
 *
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

    private CRC32 crc32 = new CRC32();

    @Override
    public ByteArrayOutputStream encode(final Message readCollection) throws IOException {
        if (readCollection == null) {
            return null;
        }
        final ByteArrayOutputStream result = new ByteArrayOutputStream();
        final DataOutputStream completeChunkData = new DataOutputStream(result);
        final ByteArrayOutputStream hybridStreamBytes = new ByteArrayOutputStream();
        final Message reducedProtoBuff = handler.compressCollection(readCollection, hybridStreamBytes);

        final int hybridStreamSize = hybridStreamBytes.size();
        final byte[] bytes = hybridStreamBytes.toByteArray();

        crc32.reset();
        crc32.update(bytes);
        final int crcChecksum = (int) crc32.getValue();
        completeChunkData.writeInt(hybridStreamSize);
        completeChunkData.writeInt(crcChecksum);
        completeChunkData.write(bytes);

        final ByteArrayOutputStream out = gzipCodec.encode(reducedProtoBuff);

        final byte[] gzipBytes = out.toByteArray();
        final int gzipBytesSize = gzipBytes.length;
        completeChunkData.write(gzipBytes);
        completeChunkData.flush();
        if (debug && chunkIndex % 100 == 0) {

            //TODO remove compression of original collection. Only useful for stat collection
            int originalGzipSize = gzipCodec.encode(readCollection).toByteArray().length;

            final int gain = originalGzipSize - (gzipBytesSize + hybridStreamSize);
            LOG.info(String.format("compressed size=%d gzip size=%d (original gzip=%d) percent compressed/(compressed+gzip) %g gain=%d, %g%% ",
                    hybridStreamSize, gzipBytesSize, originalGzipSize,
                    100d * ((double) hybridStreamSize) / (hybridStreamSize + gzipBytesSize),
                    gain, gain * 100d / originalGzipSize));

        }
        chunkIndex++;
        return result;
    }

    @Override
    public Message decode(final byte[] bytes) throws IOException {
        final DataInputStream completeChunkData = new DataInputStream(new ByteArrayInputStream(bytes));
        final int compressedSize = completeChunkData.readInt();
        final int storedChecksum = completeChunkData.readInt();

        final byte[] compressedBytes = new byte[compressedSize];
        final int read = completeChunkData.read(compressedBytes, 0, compressedSize);
        assert read == compressedSize : "read size must match recorded size.";
        crc32.reset();

        crc32.update(compressedBytes);
        final int computedChecksum = (int) crc32.getValue();
        if (computedChecksum != storedChecksum) {
            throw new InvalidChecksumException();
        }
        final int bytesLeft = bytes.length - 4 - compressedSize - 4;
        final byte[] leftOver = new byte[bytesLeft];
        // 8 is the number of bytes to encode the length of the compressed chunk, plus
        // the number of bytes to encode the checksum.
        System.arraycopy(bytes, 8 + compressedSize, leftOver, 0, bytesLeft);
        final Message reducedProtoBuff = gzipCodec.decode(leftOver);
        return handler.decompressCollection(reducedProtoBuff, compressedBytes);
    }

    @Override
    public int getSuggestedChunkSize() {
        return 30000;
    }

    @Override
    public void setHandler(final ProtobuffCollectionHandler handler) {
        this.handler = handler;
        gzipCodec.setHandler(handler);
    }

    @Override
    public boolean validate(byte firstByte, final DataInputStream input) {

        try {
            crc32.reset();

            final byte b = input.readByte();
            final byte c = input.readByte();
            final byte d = input.readByte();
            final int fullCodecContentSize = firstByte << 24 | (b& 0xFF) << 16 | (c & 0xFF)<< 8 | (d& 0xFF);
           // System.out.printf("read %X %X %X %X but found size=%X", firstByte, b, c, d,fullCodecContentSize);

            final int hybridContentSize = input.readInt();
            final int storedChecksum = input.readInt();
            if (fullCodecContentSize < 0) {
                return false;
            }

            if (hybridContentSize < 0) {
                return false;
            }

            final byte[] bytes = new byte[hybridContentSize];
            int totalRead = 0;
            int offset = 0;
            while (totalRead < hybridContentSize) {
                final int numRead = input.read(bytes, offset, hybridContentSize - totalRead);
                if (numRead == -1) {
                    break;
                }
                totalRead += numRead;
                offset += numRead;

            }
            if (totalRead != hybridContentSize) {
                return false;
            }
            crc32.update(bytes);
            final int computedChecksum = (int) crc32.getValue();
            return computedChecksum == storedChecksum;
        } catch (IOException e) {
            return false;
        }

    }


}
