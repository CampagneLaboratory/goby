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

import com.google.protobuf.ByteString;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.junit.Before;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 1/16/12
 *         Time: 10:58 AM
 */
public class TestReadCodecImpl {
    class ReadExample {
        int readIndex;
        int readLength;
        String sequence;
        String quality;
        String sequencePair;

        ReadExample(int readIndex, int readLength, String sequence, String quality, String sequencePair, String qualityPair) {
            this.readIndex = readIndex;
            this.readLength = readLength;
            this.sequence = sequence;
            this.quality = quality;
            this.sequencePair = sequencePair;
            this.qualityPair = qualityPair;
        }

        String qualityPair;
    }

    ReadExample[] examples = new ReadExample[]{
            new ReadExample(0, 10, "ACTGATTCAC", "hhhhhhhhhh", null, null),
            new ReadExample(1, 10, "ACTGATTCAC", "hhhhhhhhhh", null, null)
    };
    ObjectArrayList<Reads.ReadEntry.Builder> builtEntries;

    @Before
    public void setup() {
        builtEntries = new ObjectArrayList<Reads.ReadEntry.Builder>();
        for (ReadExample entry : examples) {
            Reads.ReadEntry.Builder readBuilder = Reads.ReadEntry.newBuilder();
            readBuilder.setReadIndex(entry.readIndex);
            readBuilder.setReadLength(entry.readLength);
            if (entry.sequence != null) {
                readBuilder.setSequence(convertSequence(entry.sequence));
            }
            if (entry.quality != null) {
                readBuilder.setQualityScores(convertQuality(entry.quality));
            }
            if (entry.sequencePair != null) {
                readBuilder.setSequencePair(convertSequence(entry.sequencePair));
            }
            if (entry.qualityPair != null) {
                readBuilder.setQualityScoresPair(convertQuality(entry.qualityPair));
            }

            builtEntries.add(readBuilder);
        }

    }

    private ByteString convertQuality(String quality) {
        byte[] result = new byte[quality.length()];

        for (int i = 0; i < quality.length(); i++) {
            result[i] = (byte) quality.charAt(i);
        }
        return ByteString.copyFrom(result);
    }

    private ByteString convertSequence(String sequence) {
        return ReadsWriterImpl.encodeSequence(sequence, new byte[sequence.length()]);
    }

    @Test
    public void testCodec() {
        ReadCodec codec = new ReadCodecImpl();
        for (Reads.ReadEntry.Builder entry : builtEntries) {

            Reads.ReadEntry.Builder compressed = codec.encode(entry);
            codec.newChunk();
            Reads.ReadEntry.Builder decompressed = codec.decode(compressed.build());
            assertEquals(entry.build().toString(), decompressed.build().toString());
        }
    }
}
