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

import edu.cornell.med.icb.goby.alignments.AlignmentCollectionHandler;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.junit.Before;
import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;

/**
 * @author Fabien Campagne
 *         Date: 3/3/12
 *         Time: 3:40 PM
 */
public class TestAlignmentChunkCodec1 {
    @Test
    public void nullCollection() throws IOException {
        final AlignmentChunkCodec1 codec = new AlignmentChunkCodec1();
        assertEquals(null, codec.encode(null));
    }

    @Test
    public void smallCollection() throws IOException {
        final AlignmentChunkCodec1 codec = new AlignmentChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        final AlignmentReader reader = new AlignmentReaderImpl("test-data/alignments/mantis-1355/mantis-1355.entries");
        final Alignments.AlignmentCollection.Builder collectionBuilder = Alignments.AlignmentCollection.newBuilder();
        for (Alignments.AlignmentEntry entry : reader) {
            collectionBuilder.addAlignmentEntries(entry);
        }

        final ByteArrayOutputStream encoded = codec.encode(collectionBuilder.build());
        assertNotNull(encoded);
        // assertEquals(collectionBuilder.getAlignmentEntriesCount(), encoded.);
    }


    // @Test
    // TODO enable to debug decoding:
    public void testRoundTrip() throws IOException {
        final AlignmentChunkCodec1 codec = new AlignmentChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = buildCollection();

        final ByteArrayOutputStream encoded = codec.encode(collection.build());
        Alignments.AlignmentCollection decodedCollection = (Alignments.AlignmentCollection) codec.decode(encoded.toByteArray());

        for (int i=0;i<decodedCollection.getAlignmentEntriesCount(); i++) {

            Alignments.AlignmentEntry decodedEntry = decodedCollection.getAlignmentEntries(i);
            final String expected = builtEntries.get(i).build().toString();
            final String actual = decodedEntry.toString();
            assertEquals("entry " + i, expected, actual);
        }
    }

    private Alignments.AlignmentCollection.Builder buildCollection() {
        final Alignments.AlignmentCollection.Builder collectionBuilder = Alignments.AlignmentCollection.newBuilder();
        for (Alignments.AlignmentEntry.Builder entry : builtEntries) {
            collectionBuilder.addAlignmentEntries(entry);
        }
        return collectionBuilder;
    }

    @Test
    public void testMore() throws IOException {
        final AlignmentChunkCodec1 codec = new AlignmentChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        final AlignmentReader reader = new AlignmentReaderImpl("test-data/seq-var-test/kevin-synth/sorted-seq-var-reads-gsnap.entries");
        final Alignments.AlignmentCollection.Builder collectionBuilder = Alignments.AlignmentCollection.newBuilder();
        for (Alignments.AlignmentEntry entry : reader) {
            collectionBuilder.addAlignmentEntries(entry);
        }

        final ByteArrayOutputStream encoded = codec.encode(collectionBuilder.build());
        assertNotNull(encoded);
        // assertEquals(collectionBuilder.getAlignmentEntriesCount(), encoded.);
    }
    // @Test
    public void testLarge() throws IOException {
        final AlignmentChunkCodec1 codec = new AlignmentChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        final AlignmentReader reader = new AlignmentReaderImpl("/data/rrbs/EMNWFIL.entries");
        final Alignments.AlignmentCollection.Builder collectionBuilder = Alignments.AlignmentCollection.newBuilder();
        int i = 0;
        for (Alignments.AlignmentEntry entry : reader) {
            collectionBuilder.addAlignmentEntries(entry);
            if (i++ > 10000) {
                break;
            }
        }

        final ByteArrayOutputStream encoded = codec.encode(collectionBuilder.build());
        assertNotNull(encoded);
        // assertEquals(collectionBuilder.getAlignmentEntriesCount(), encoded.);
    }


    class AlignmentExample {
        int position;
        int mappingQuality;
        int query_index;

        String var1_to;
        String var1_from;
        int var1_position;
        int var1_readIndex;

        String var2_to;
        String var2_from;
        int var2_position;
        int var2_readIndex;


        AlignmentExample(int position, int mappingQuality, int query_index, String var1_to, String var1_from, int var1_position, int var1_readIndex,
                         String var2_to, String var2_from, int var2_position, int var2_readIndex) {

            this.position = position;
            this.mappingQuality = mappingQuality;

            this.var1_to = var1_to;
            this.var1_from = var1_from;
            this.var1_position = var1_position;
            this.var1_readIndex = var1_readIndex;

            this.var2_to = var2_to;
            this.var2_from = var2_from;
            this.var2_position = var2_position;
            this.var2_readIndex = var2_readIndex;

            //required fields that are uncompressed:
            this.query_index = query_index;
        }
    }

    AlignmentExample[] examples = new AlignmentExample[]{
            new AlignmentExample(1, 33, 1, "TA", "GG", 27, 31, "--", "TA", 3, 1000),
            new AlignmentExample(1, 9, 1, "..", "AA", 42, 24, "T-", "G.", 5, 1001)
    };
    ObjectArrayList<Alignments.AlignmentEntry.Builder> builtEntries;

    @Before
    public void setup() {
        builtEntries = new ObjectArrayList<Alignments.AlignmentEntry.Builder>();
        for (AlignmentExample entry : examples) {
            Alignments.AlignmentEntry.Builder alignmentBuilder = Alignments.AlignmentEntry.newBuilder();
            alignmentBuilder.setPosition(entry.position);
            alignmentBuilder.setMappingQuality(entry.mappingQuality);
            alignmentBuilder.setQueryIndex(entry.query_index);

            Alignments.SequenceVariation.Builder sequenceVariation1 = Alignments.SequenceVariation.newBuilder();
            sequenceVariation1.setFrom(entry.var1_from);
            sequenceVariation1.setTo(entry.var1_to);
            sequenceVariation1.setPosition(entry.var1_position);
            sequenceVariation1.setReadIndex(entry.var1_readIndex);
            alignmentBuilder.addSequenceVariations(sequenceVariation1.build());

            Alignments.SequenceVariation.Builder sequenceVariation2 = Alignments.SequenceVariation.newBuilder();
            sequenceVariation2.setFrom(entry.var2_from);
            sequenceVariation2.setTo(entry.var2_to);
            sequenceVariation2.setPosition(entry.var2_position);
            sequenceVariation2.setReadIndex(entry.var2_readIndex);
            alignmentBuilder.addSequenceVariations(sequenceVariation2.build());

            builtEntries.add(alignmentBuilder);
        }

    }


}
