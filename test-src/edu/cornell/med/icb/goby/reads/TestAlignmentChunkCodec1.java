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
import edu.cornell.med.icb.goby.alignments.perms.QueryIndexPermutation;
import edu.cornell.med.icb.goby.alignments.perms.QueryIndexPermutationInterface;
import edu.cornell.med.icb.goby.compression.HybridChunkCodec1;
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
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        assertEquals(null, codec.encode(null));
    }

    @Test
    public void smallCollection() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
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


    @Test

    public void testRoundTrip() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = buildCollection(examples);


        testRoundTripWithBuiltEntries(codec, collection);
    }

    AlignmentExample[] examplesWithDuplicates = {
            new AlignmentExample(0, 1, 33, 1, "TA", "GG", 27, 31, "--", "TA", 3, 1000),    // R1
            new AlignmentExample(0, 2, 33, 2, "TA", "GG", 27, 31, "--", "TA", 3, 1000),    // R1
            new AlignmentExample(0, 3, 33, 3, "TA", "GG", 27, 31, "--", "TA", 3, 1000),    // R1
            new AlignmentExample(0, 15, 33, 4, "TA", "GG", 27, 31, "--", "TA", 3, 1000),   // R1
            new AlignmentExample(0, 101, 9, 5, "..", "AA", 42, 24, "T-", "G.", 5, 1001),     // R2
            new AlignmentExample(0, 112, 9, 6, "..", "AA", 42, 24, "T-", "G.", 5, 1001),    // R2
            new AlignmentExample(0, 112, 9, 7, "..", "AA", 42, 24, "T-", "G.", 5, 1001),    // R2
            new AlignmentExample(1, 1111, 9, 8, "AA--AA", "CCCC", 13, 24, "T-", "G.", 5, 3)   //R3
    };

    @Test
    public void testRoundTripMultipleEntriesDuplicate() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = buildCollection(examplesWithDuplicates);
        //  System.out.println(collection.build().toString());
        assertRoundTripMatchExpected(codec, collection);
    }

    @Test
    public void roundTripMore() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = loadCollection("test-data/seq-var-test/kevin-synth/sorted-seq-var-reads-gsnap.entries");

        assertRoundTripMatchExpected(codec, collection);
    }

    //@Test
    // will not run on server.
    public void roundTripLarge() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = loadCollection("/data/rrbs/EMNWFIL.entries", 0, 100);

        assertRoundTripMatchExpected(codec, collection, false);
    }

    @Test

    public void roundTripExamplePairedEnd() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        final Alignments.AlignmentCollection.Builder collection = loadCollection("test-data/bam/Example.entries",0,1000);

        assertRoundTripMatchExpected(codec, collection);
    }

    // @Test
    // will not run on server.
    public void roundTripPairedEnd() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = loadCollection("/data/CRAM/VJDQTEI-C1.entries", 600, 900);

        assertRoundTripMatchExpected(codec, collection);
    }

    //  @Test
    // will not run on server.
    public void roundTripBug() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = loadCollection("/data/CRAM/WZLFUIH-paper-combined-NA18853.entries", 10000 * 117 - 5000, 10000);

        assertRoundTripMatchExpected(codec, collection);
    }

    private Alignments.AlignmentCollection.Builder loadCollection(String filename) throws IOException {
        return loadCollection(filename, 0, Integer.MAX_VALUE);
    }

    private Alignments.AlignmentCollection.Builder loadCollection(String filename, int firstElementToLoad, int maxElementsToLoad) throws IOException {
        final Alignments.AlignmentCollection.Builder collectionBuilder = Alignments.AlignmentCollection.newBuilder();
        AlignmentReaderImpl reader = new AlignmentReaderImpl(filename);
        QueryIndexPermutationInterface permutator = new QueryIndexPermutation(filename);
        try {
            int counter = 0;
            for (Alignments.AlignmentEntry entry : reader) {
                if (counter >= firstElementToLoad) {
                    collectionBuilder.addAlignmentEntries(permutator.makeSmallIndices(entry));
                    if (counter > maxElementsToLoad) {
                        break;
                    }
                }
                counter++;
            }
            return collectionBuilder;
        } finally {
            reader.close();
        }
    }

    private void testRoundTripWithBuiltEntries(HybridChunkCodec1 codec, Alignments.AlignmentCollection.Builder collection) throws IOException {
        final ByteArrayOutputStream encoded = codec.encode(collection.build());
        Alignments.AlignmentCollection decodedCollection = (Alignments.AlignmentCollection) codec.decode(encoded.toByteArray());
        Alignments.AlignmentCollection.Builder expected = Alignments.AlignmentCollection.newBuilder();
        for (final Alignments.AlignmentEntry.Builder entryBuilder : builtEntries) {
            expected.addAlignmentEntries(entryBuilder);
        }
        assertEquals("collection", expected.build().toString(), decodedCollection.toString());

    }

    private void assertRoundTripMatchExpected(HybridChunkCodec1 codec, Alignments.AlignmentCollection.Builder expected) throws IOException {
        assertRoundTripMatchExpected(codec, expected, true);
    }

    private void assertRoundTripMatchExpected(HybridChunkCodec1 codec, Alignments.AlignmentCollection.Builder expected, boolean doAssert) throws IOException {
        final ByteArrayOutputStream encoded = codec.encode(expected.build());
        Alignments.AlignmentCollection decodedCollection = (Alignments.AlignmentCollection) codec.decode(encoded.toByteArray());

        if (doAssert) {
            assertEquals("collection", expected.build().toString(), decodedCollection.toString());
        }

    }

    private Alignments.AlignmentCollection.Builder buildCollection(AlignmentExample[] builtEntries) {
        final Alignments.AlignmentCollection.Builder collectionBuilder = Alignments.AlignmentCollection.newBuilder();
        for (Alignments.AlignmentEntry.Builder entry : buildEntriesCollection(builtEntries)) {
            collectionBuilder.addAlignmentEntries(entry);
        }
        return collectionBuilder;
    }


    @Test
    public void testMore() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
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
    // will not run on server.
    public void testLarge() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
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
        int targetIndex;

        String var1_to;
        String var1_from;
        int var1_position;
        int var1_readIndex;

        String var2_to;
        String var2_from;
        int var2_position;
        int var2_readIndex;


        AlignmentExample(int targetIndex, int position, int mappingQuality, int query_index, String var1_to, String var1_from, int var1_position, int var1_readIndex,
                         String var2_to, String var2_from, int var2_position, int var2_readIndex) {

            this.targetIndex = targetIndex;
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
            new AlignmentExample(0, 1, 33, 1, "TA", "GG", 27, 31, "--", "TA", 3, 1000),
            new AlignmentExample(0, 1, 9, 1, "..", "AA", 42, 24, "T-", "G.", 5, 1001),
            new AlignmentExample(1, 1, 9, 1, "AA--AA", "CCCC", 13, 24, "T-", "G.", 5, 3)
    };
    ObjectArrayList<Alignments.AlignmentEntry.Builder> builtEntries;

    @Before
    public void setup() {
        builtEntries = buildEntriesCollection(examples);


    }

    private ObjectArrayList<Alignments.AlignmentEntry.Builder> buildEntriesCollection(AlignmentExample[] examples) {
        ObjectArrayList<Alignments.AlignmentEntry.Builder> list = new ObjectArrayList<Alignments.AlignmentEntry.Builder>();
        for (AlignmentExample entry : examples) {
            Alignments.AlignmentEntry.Builder alignmentBuilder = Alignments.AlignmentEntry.newBuilder();
            alignmentBuilder.setPosition(entry.position);
            alignmentBuilder.setMappingQuality(entry.mappingQuality);
            alignmentBuilder.setQueryIndex(entry.query_index);
            alignmentBuilder.setMatchingReverseStrand(true);
            alignmentBuilder.setTargetIndex(entry.targetIndex);
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

            list.add(alignmentBuilder);
        }
        return list;
    }


}
