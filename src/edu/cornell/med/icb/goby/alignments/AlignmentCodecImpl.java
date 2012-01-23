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

package edu.cornell.med.icb.goby.alignments;

import com.google.protobuf.ByteString;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.io.FastByteArrayOutputStream;
import it.unimi.dsi.io.InputBitStream;
import it.unimi.dsi.io.OutputBitStream;

import java.io.IOException;

/**
 * A codec to compress alignments. Alignments are written to the special
 * compressed_data protocol buffer field. This class provides a proof of
 * principle for combining data-specific compression schemes with the
 * flexibility of protocol buffer.
 *
 * @author Michael Meyer
 *         Date: 1/5/12
 *         Time: 2:30 PM
 */
public class AlignmentCodecImpl implements AlignmentCodec {
    public static final int CODEC_REGISTRATION_CODE = 2;
    private int Encode_previousPosition;
    private int Decode_previousPosition;
    private int Encode_previousMappingQuality;
    private int Decode_previousMappingQuality;

    @Override
    public Alignments.AlignmentEntry.Builder encode(final Alignments.AlignmentEntry.Builder source) {
        final Alignments.AlignmentEntry.Builder result = Alignments.AlignmentEntry.newBuilder();

        FastByteArrayOutputStream os = new FastByteArrayOutputStream();
        // compress sequence variation of the alignment:

//        OutputBitStream out = new DebugOutputBitStream(new OutputBitStream(os));
        OutputBitStream out = new OutputBitStream(os);
        result.mergeFrom(source.build());
        try {
//             write the codec registration code first as one byte:
            out.writeInt(CODEC_REGISTRATION_CODE, 8);

            out.writeUnary(source.getSequenceVariationsCount());
            if (source.getSequenceVariationsCount() != 0) {
                for (int i = 0; i < source.getSequenceVariationsCount(); i++) {
                    compressSequenceVariation(source.getSequenceVariations(i), out);
                }
            }

            if (source.hasPosition()) {
                out.writeBit(true);
                compressDelta(source.getPosition(), Encode_previousPosition, "P", out);
                Encode_previousPosition = source.getPosition();
            } else {
                out.writeBit(false);
            }

            if (source.hasMappingQuality()) {
                out.writeBit(true);
                compressDelta(source.getMappingQuality(), Encode_previousMappingQuality, "MQ", out);
                Encode_previousMappingQuality = source.getMappingQuality();


            } else {
                out.writeBit(false);
            }

            out.flush();
            final ByteString compressedData = ByteString.copyFrom(os.array, 0, (int) os.length());

            result.setCompressedData(compressedData);
            result.clearSequenceVariations();
            result.clearPosition();
            result.clearMappingQuality();
            return result;

        } catch (IOException e) {
            // an error occurred compressing the alignment. Return null to indicate that we cannot handle the input.
            return null;
        }

    }


    private void writeBit(OutputBitStream out, boolean bit) throws IOException {
        out.writeBit(bit);
    }

    private void writeUnary(OutputBitStream out, int num) throws IOException {
        out.writeUnary(num);
    }


    private void compressDelta(final int current, final int previous, String source, final OutputBitStream out) throws IOException {
        //for compressing sequences of numbers with small variations, like position, index
        int difference = current - previous;

        if (difference >= 0) {
            //delta coding is incapable of encoding zero, (when decoding, subtract 1 from all positive numbers)
            out.writeBit(true);
            out.writeDelta(difference + 1);
        } else {
            out.writeBit(false);
            out.writeDelta(Math.abs(difference));
        }

    }

    private void compressSequenceVariation(final edu.cornell.med.icb.goby.alignments.Alignments.SequenceVariation sequenceVariation, final OutputBitStream out) throws IOException {
        int[] bits;
        if (true) {
            String from = sequenceVariation.getFrom();
            //encode length of from field using Unary coding - usually 1 bit
            writeUnary(out, from.length());
            for (int s = 0; s < from.length(); s++) {
                bits = codeBase(from.charAt(s));
                for (int r = 0; r < bits.length; r++) {
                    out.writeBit(bits[r]);
                }
            }
        }

        if (true) {
            String to = sequenceVariation.getTo();
            writeUnary(out, to.length());
            for (int s = 0; s < to.length(); s++) {
                bits = codeBase(to.charAt(s));
                for (int r = 0; r < bits.length; r++) {
                    out.writeBit(bits[r]);
                }
            }
        }

        if (true) {
            out.writeDelta(sequenceVariation.getPosition() + 1);
        }

        if (true) {
            out.writeDelta(sequenceVariation.getReadIndex() + 1);
        }

        if (sequenceVariation.hasToQuality()) {
            ByteString toQuality = sequenceVariation.getToQuality();
            out.writeUnary(toQuality.size());
            for (int j = 0; j < toQuality.size(); j++) {
                out.writeDelta(((int) toQuality.byteAt(j)) + 1);
            }
        } else {
            out.writeUnary(0);
        }

    }


    private int[] codeBase(final char base) {
        switch (base) {
            case 'A':
                return new int[]{0, 0, 0};
            case 'C':
                return new int[]{0, 0, 1};
            case 'T':
                return new int[]{0, 1, 0};
            case 'G':
                return new int[]{0, 1, 1};
            case 'N':
                return new int[]{1, 0, 0};
            case '-':
                return new int[]{1, 0, 1};
            case '.':
                return new int[]{1, 1, 0};
        }
        return new int[]{1, 1, 1};
    }

    private char decodeBase(int[] decode) {
        switch (decode[0]) {
            case 0:
                switch (decode[1]) {
                    case 0:
                        switch (decode[2]) {
                            case 0:
                                return 'A';
                            case 1:
                                return 'C';
                        }
                    case 1:
                        switch (decode[2]) {
                            case 0:
                                return 'T';
                            case 1:
                                return 'G';
                        }
                }
            case 1:
                switch (decode[1]) {
                    case 0:
                        switch (decode[2]) {
                            case 0:
                                return 'N';
                            case 1:
                                return '-';
                        }
                    case 1:
                        switch (decode[2]) {
                            case 0:
                                return '.';
                            case 1:
                                return '9';
                        }
                }
        }
        return '9';
    }

    @Override
    public Alignments.AlignmentEntry.Builder decode(Alignments.AlignmentEntry source) {
        if (!source.hasCompressedData()) {
            return null;
        }
        InputBitStream input = new InputBitStream(source.getCompressedData().toByteArray());
        try {
            final int codecRegistrationStored = input.readInt(8);
            if (codecRegistrationStored != CODEC_REGISTRATION_CODE) {
                // this alignment is compressed by a different codec, indicate that we cannot handle it.
                return null;
            }
            final Alignments.AlignmentEntry.Builder result = Alignments.AlignmentEntry.newBuilder();
            // get any other fields the codec does not handle:
            result.mergeFrom(source);

            int sequenceVariationCount = input.readUnary();
            if (sequenceVariationCount > 0) {

                for (int i = 0; i < sequenceVariationCount; i++) {
                    result.addSequenceVariations(decodeSequenceVariations(input));
                }
            }

            //Decode Position (1 means position was encoded)
            int curPosition;
            if (input.readBit() == 1) {

                curPosition = decodeDelta(Decode_previousPosition, input);
                result.setPosition(curPosition);
                Decode_previousPosition = curPosition;

            }

            int curMappingQuality;
            if (input.readBit() == 1) {

                curMappingQuality = decodeDelta(Decode_previousMappingQuality, input);
                result.setMappingQuality(curMappingQuality);
                Decode_previousMappingQuality = curMappingQuality;

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

        return "alignment-codec-2";
    }

    @Override
    public final void newChunk() {
        //new chunk - initialize encoders/decoders
        Encode_previousPosition = 0;
        Decode_previousPosition = 0;
        Encode_previousMappingQuality = 0;
        Decode_previousMappingQuality = 0;
    }

    public AlignmentCodecImpl() {
        newChunk();
    }

    private int decodeDelta(int previous, InputBitStream input) throws IOException {

        int result;
        //find if the difference was positive or negative (pos=1, neg=0)
        if (input.readBit() == 1) {
            //positive numbers have been incremented by 1, so subtract one from the difference to regain the number
            result = previous + input.readDelta() - 1;
        } else {
            result = previous - input.readDelta();
        }

        return result;

    }

    private edu.cornell.med.icb.goby.alignments.Alignments.SequenceVariation decodeSequenceVariations(InputBitStream input) throws IOException {

        final Alignments.SequenceVariation.Builder currentSequenceVariation = Alignments.SequenceVariation.newBuilder();

        //Decode <from>

        int fromLength = input.readUnary();

        if (fromLength > 0) {
            StringBuffer from = new StringBuffer();
            for (int j = 0; j < fromLength; j++) {
                int bit1 = input.readBit();
                int bit2 = input.readBit();
                int bit3 = input.readBit();
                from.append(decodeBase(new int[]{bit1, bit2, bit3}));
            }
            currentSequenceVariation.setFrom(from.toString());
        }

        //Decode <to>

        int toLength = input.readUnary();

        if (toLength > 0) {
            StringBuffer to = new StringBuffer();
            for (int j = 0; j < toLength; j++) {
                int bit1 = input.readBit();
                int bit2 = input.readBit();
                int bit3 = input.readBit();
                to.append(decodeBase(new int[]{bit1, bit2, bit3}));
            }
            currentSequenceVariation.setTo(to.toString());
        }

        //Decode <position>

        currentSequenceVariation.setPosition(input.readDelta() - 1);

        //Decode <read-index>

        currentSequenceVariation.setReadIndex(input.readDelta() - 1);

        //Decode <to-quality>?

        int toQualityLength = input.readUnary();
        if (toQualityLength > 0) {
            ByteArrayList buffer = new ByteArrayList();
            for (int i = 0; i > toQualityLength; i++) {
                buffer.add((byte) (input.readDelta() - 1));
            }
            currentSequenceVariation.setToQuality(ByteString.copyFrom(buffer.toByteArray()));

        }

        return currentSequenceVariation.build();


    }

    @Override
    public byte registrationCode() {
        return CODEC_REGISTRATION_CODE;
    }
}
    
