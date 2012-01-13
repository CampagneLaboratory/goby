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
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.io.FastByteArrayOutputStream;
import it.unimi.dsi.io.InputBitStream;
import it.unimi.dsi.io.OutputBitStream;
import it.unimi.dsi.mg4j.io.ArithmeticCoder;
import it.unimi.dsi.mg4j.io.ArithmeticDecoder;

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
    private ArithmeticCoder sequenceCoder;
    private ArithmeticCoder qualityScoreCoder;
    public static final int CODEC_REGISTRATION_CODE = 1;
    FastByteArrayOutputStream os = new FastByteArrayOutputStream();
             OutputBitStream out = new OutputBitStream(os);

    @Override
    public Reads.ReadEntry.Builder encode(final Reads.ReadEntry.Builder source) {
        final Reads.ReadEntry.Builder result = Reads.ReadEntry.newBuilder();

        // compress sequence of the read:

        result.mergeFrom(source.build());
        try {
            out.flush();
            os.reset();
            // write the codec registration code first as one byte:
            out.writeInt(CODEC_REGISTRATION_CODE, 8);

            if (source.hasSequence()) {
                writeBit(out, true);
                compressSequence(source.getSequence(), out);
                result.clearSequence();
            } else {
                writeBit(out, false);
            }

            if (source.hasQualityScores()) {
                writeBit(out, true);
                compressQuality(source.getQualityScores(), out);
                result.clearQualityScores();
            } else {
                writeBit(out, false);
            }

            if (source.hasSequencePair()) {
                writeBit(out, true);
                compressSequence(source.getSequencePair(), out);
                result.clearSequencePair();
            } else {
                writeBit(out, false);
            }
            if (source.hasQualityScoresPair()) {
                writeBit(out, true);
                compressQuality(source.getQualityScoresPair(), out);
                result.clearQualityScoresPair();
            } else {
                writeBit(out, false);
            }
            out.flush();
            final ByteString compressedData = ByteString.copyFrom(os.array, 0, (int)os.length());

            result.setCompressedData(compressedData);

            return result;
        } catch (IOException e) {
            // an error occurred compressing the read. Return null to indicate that we cannot handle the input.
            return null;
        }

    }

    private void writeBit(OutputBitStream out, boolean bit) throws IOException {
        out.writeBit(bit);
    }

    private void compressQuality(final ByteString qualityScores, final OutputBitStream out) throws IOException {
        for (int i = 0; i < qualityScores.size(); i++) {
            qualityScoreCoder.encode(qualityScores.byteAt(i), out);
        }
        qualityScoreCoder.flush(out);
    }

    private void compressSequence(final ByteString sequence, final OutputBitStream out) throws IOException {
        for (int i = 0; i < sequence.size(); i++) {
            sequenceCoder.encode(codeBase(sequence.byteAt(i)), out);
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
    public Reads.ReadEntry.Builder decode(Reads.ReadEntry source) {
        if (!source.hasCompressedData()) {
            return null;
        }
        InputBitStream input = new InputBitStream(source.getCompressedData().toByteArray());
        try {
            final int codecRegistrationStored = input.readInt(8);
            if (codecRegistrationStored != CODEC_REGISTRATION_CODE) {
                // this read is compressed by a different codec, indicate that we cannot handle it.
                return null;
            }
            final Reads.ReadEntry.Builder result = Reads.ReadEntry.newBuilder();
            // get any other fields the codec does not handle:
            result.mergeFrom(source);

            if (input.readBit() == 1) {
                // sequence was stored, decode it.
                final ByteString sequence = decodeSequence(input, source.getReadLength());
                result.setSequence(sequence);
            }
            if (input.readBit() == 1) {
                // quality scores were stored, decode it.
                final ByteString qual = decodeQualityScore(input, source.getReadLength());
                result.setQualityScores(qual);
            }
            if (input.readBit() == 1) {
                // sequence pair was stored, decode it.
                final ByteString sequencePair = decodeSequence(input, source.getReadLength());
                result.setSequencePair(sequencePair);
            }
            if (input.readBit() == 1) {
                // quality score pair was stored, decode it.
                final ByteString qualPair = decodeSequence(input, source.getReadLength());
                result.setQualityScoresPair(qualPair);
            }
            // the compressed data was decoded, remove it:
            result.clearCompressedData();
            return result;

        } catch (IOException e) {
            // An exception occurred decoding compressed data, return null to indicate that we cannot handle this read.
            return null;
        }
    }

    @Override
    public String name() {

        return "read-codec-1";
    }

    @Override
    public byte registrationCode() {
        return CODEC_REGISTRATION_CODE;
    }

    @Override
    public final void newChunk() {
        sequenceCoder = new ArithmeticCoder(5);
        sequenceDecoder = new ArithmeticDecoder(5);
        qualityScoreCoder = new ArithmeticCoder(255);
        qualityScoreDecoder = new ArithmeticDecoder(255);
    }

    public ReadCodecImpl() {
        newChunk();
    }

    private ArithmeticDecoder sequenceDecoder;
    private ArithmeticDecoder qualityScoreDecoder;

    private ByteString decodeSequence(InputBitStream input, int readLength) throws IOException {
        ByteArrayList buffer = new ByteArrayList(readLength);
        for (int i = 0; i < readLength; i++) {
            buffer.add(decodeBase(sequenceDecoder.decode(input)));
        }
        return ByteString.copyFrom(buffer.toByteArray());
    }

    private ByteString decodeQualityScore(InputBitStream input, int readLength) throws IOException {
        ByteArrayList buffer = new ByteArrayList(readLength);
        for (int i = 0; i < readLength; i++) {
            buffer.add((byte) qualityScoreDecoder.decode(input));
        }
        return ByteString.copyFrom(buffer.toByteArray());
    }

}
