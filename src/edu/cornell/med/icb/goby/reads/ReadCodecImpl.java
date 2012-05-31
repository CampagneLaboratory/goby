/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

import com.google.protobuf.ByteString;
import edu.cornell.med.icb.goby.algorithmic.compression.FastArithmeticCoder;
import edu.cornell.med.icb.goby.algorithmic.compression.FastArithmeticDecoder;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.io.FastByteArrayOutputStream;
import it.unimi.dsi.io.InputBitStream;
import it.unimi.dsi.io.OutputBitStream;

import java.io.IOException;

/**
 * A codec to compress sequences and quality scores. Sequences and quality scores are written to the special
 * compressed_data protocol buffer field. This class provides a proof of principle for combining data-specific
 * compression schemes with the flexibility of protocol buffer.
 *
 * @author Fabien Campagne
 *         Date: 8/19/11
 *         Time: 2:00 PM
 */
public class ReadCodecImpl implements ReadCodec {
    private FastArithmeticCoder sequenceCoder;
    private FastArithmeticCoder qualityScoreCoder;
    public static final int CODEC_REGISTRATION_CODE = 1;
    private boolean isFirstCode=true;
    private boolean isFirstDecode;

    @Override
    public String name() {

        return "read-codec-1";
    }

    @Override
    public byte registrationCode() {
        return CODEC_REGISTRATION_CODE;
    }

    int written = 0;

    @Override
    public Reads.ReadEntry.Builder encode(final Reads.ReadEntry.Builder source) {
        final Reads.ReadEntry.Builder result = Reads.ReadEntry.newBuilder();

        // compress sequence of the read:

        result.mergeFrom(source.build());
        try {
            FastByteArrayOutputStream os = new FastByteArrayOutputStream();
            OutputBitStream out = /*new DebugOutputBitStream(*/new OutputBitStream(os)/*)*/;

            if (isFirstCode) {
                // write the codec registration code first as one byte:
                out.writeInt(CODEC_REGISTRATION_CODE, 8);
                isFirstCode=false;
            }
            final int readLength = source.getReadLength();

            if (source.hasSequence()) {

                writeBit(out, true);
                compressSequence(source.getSequence(), out, readLength);
                result.clearSequence();
            } else {
                writeBit(out, false);
            }

            if (source.hasQualityScores()) {
                writeBit(out, true);
                compressQuality(source.getQualityScores(), out, readLength);
                result.clearQualityScores();
            } else {
                writeBit(out, false);
            }

            if (source.hasSequencePair()) {
                writeBit(out, true);
                compressSequence(source.getSequencePair(), out, readLength);
                result.clearSequencePair();
            } else {
                writeBit(out, false);
            }
            if (source.hasQualityScoresPair()) {

                writeBit(out, true);
                compressQuality(source.getQualityScoresPair(), out, readLength);
                result.clearQualityScoresPair();
            } else {
                writeBit(out, false);
            }
            out.close();
            final ByteString compressedData = ByteString.copyFrom(os.array, 0, (int) os.length());

            result.setCompressedData(compressedData);

            return result;
        } catch (IOException e) {
            // an error occurred compressing the read. Return null to indicate that we cannot handle the input.
            return null;
        }

    }

    private void writeBit(OutputBitStream out, boolean bit) throws IOException {
        final int i = out.writeBit(bit);
        assert i == 1;
    }

    private void compressQuality(final ByteString qualityScores, final OutputBitStream out, int readLength) throws IOException {
        qualityScoreCoder.reset();
        for (int i = 0; i < readLength; i++) {
            final byte x = qualityScores.byteAt(i);

            qualityScoreCoder.encode(x, out);
        }
        qualityScoreCoder.flush(out);

    }

    private void compressSequence(final ByteString sequence, final OutputBitStream out, int readLength) throws IOException {
        sequenceCoder.reset();
        for (int i = 0; i < readLength; i++) {
            final int x = codeBase(sequence.byteAt(i));

            sequenceCoder.encode(x, out);
        }
        sequenceCoder.flush(out);

    }

    private int codeBase(final byte base) {
        switch (base) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'T':
                return 2;
            case 'G':
                return 3;
            case 'N':
                return 4;
        }
        return -1;
    }

    private byte decodeBase(int decode) {
        switch (decode) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'T';
            case 3:
                return 'G';
            case 4:
                return 'N';
        }
        return -1;
    }

    @Override
    public Reads.ReadEntry.Builder decode(final Reads.ReadEntry source) {
        if (!source.hasCompressedData()) {
            return null;
        }
        byte[] bytes = new byte[source.getCompressedData().size() + 40];
        source.getCompressedData().copyTo(bytes, 0);
        final InputBitStream input =/* new DebugInputBitStream(*/new InputBitStream(bytes)/*)*/;
        try {

            if (isFirstDecode) {
                final int codecRegistrationStored = input.readInt(8);
                if (codecRegistrationStored != CODEC_REGISTRATION_CODE) {
                    // this read is compressed by a different codec, indicate that we cannot handle it.
                    return null;
                }
            }
            final Reads.ReadEntry.Builder result = Reads.ReadEntry.newBuilder();
            // get any other fields the codec does not handle:
            result.mergeFrom(source);

            final int readLength = source.getReadLength();

            debug("readLength=" + readLength);
            debug("readBits= " + input.readBits());
            if (input.readBit() == 1) {
                debug("hasSequence ");
                // sequence was stored, decode it.
                debug("readBits= " + input.readBits());
                final ByteString sequence = decodeSequence(input, readLength);
                debug("readBits= " + input.readBits());
                result.setSequence(sequence);
            }
            if (input.readBit() == 1) {
                debug("hasQual ");
                debug("readBits= " + input.readBits());
                // quality scores were stored, decode it.
                final ByteString qual = decodeQualityScore(input, readLength);
                debug("readBits= " + input.readBits());
                result.setQualityScores(qual);
            }
            if (input.readBit() == 1) {
                debug("hasSequencePair ");
                // sequence pair was stored, decode it.
                final ByteString sequencePair = decodeSequence(input, readLength);
                debug("readBits= " + input.readBits());
                result.setSequencePair(sequencePair);
            }
            if (input.readBit() == 1) {
                debug("hasQualPair ");
                // quality score pair was stored, decode it.
                final ByteString qualPair = decodeQualityScore(input, readLength);
                debug("readBits= " + input.readBits());
                result.setQualityScoresPair(qualPair);
            }
            // the compressed data was decoded, remove it:
            result.clearCompressedData();
            isFirstDecode = false;
            input.close();

            return result;

        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
            // An exception occurred decoding compressed data, return null to indicate that we cannot handle this read.
            //  return null;
        }
    }

    private void debug(final String text) {
        if (false) {
            System.out.println(text);
        }
    }


    @Override
    public final void newChunk() {
        reset();
        isFirstCode = true;
        isFirstDecode = true;
    }

    private void reset() {
        sequenceCoder = new FastArithmeticCoder(5);
        sequenceDecoder = new FastArithmeticDecoder(5);
        qualityScoreCoder = new FastArithmeticCoder(255);
        qualityScoreDecoder = new FastArithmeticDecoder(255);
    }

    public ReadCodecImpl() {
        newChunk();
    }

    private FastArithmeticDecoder sequenceDecoder;
    private FastArithmeticDecoder qualityScoreDecoder;

    ByteArrayList buffer = new ByteArrayList();

    private ByteString decodeSequence(InputBitStream input, int readLength) throws IOException {
        buffer.clear();
        sequenceDecoder.reset();
        for (int i = 0; i < readLength; i++) {
            final int decoded = sequenceDecoder.decode(input);
            buffer.add(decodeBase(decoded));
            //     System.out.printf("i=%d readBits= %d%n", i, input.readBits());
        }

        //        System.out.printf("readBits= %d%n",  input.readBits());
        sequenceDecoder.reposition(input);

        //   long flushBits= sequenceDecoder.getWindow();

        return ByteString.copyFrom(buffer.toByteArray(), 0, readLength);
    }

    private ByteString decodeQualityScore(InputBitStream input, int readLength) throws IOException {
        buffer.clear();
        qualityScoreDecoder.reset();
        //ByteArrayList buffer = new ByteArrayList(readLength);
        for (int i = 0; i < readLength; i++) {
            buffer.add((byte) qualityScoreDecoder.decode(input));
            //   System.out.printf("i=%d readBits= %d%n", i, input.readBits());
        }
        qualityScoreDecoder.reposition(input);
        return ByteString.copyFrom(buffer.toByteArray(), 0, readLength);
    }

}
