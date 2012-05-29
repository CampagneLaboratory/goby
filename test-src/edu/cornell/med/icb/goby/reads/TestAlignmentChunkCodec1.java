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
import edu.cornell.med.icb.goby.alignments.AlignmentCollectionHandler;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.perms.QueryIndexPermutation;
import edu.cornell.med.icb.goby.alignments.perms.QueryIndexPermutationInterface;
import edu.cornell.med.icb.goby.compression.HybridChunkCodec1;
import edu.cornell.med.icb.goby.compression.HybridChunkCodec2;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.junit.Before;
import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.UnsupportedEncodingException;

import static junit.framework.Assert.*;

/**
 * @author Fabien Campagne
 *         Date: 3/3/12
 *         Time: 3:40 PM
 */
public class TestAlignmentChunkCodec1 {

    @Test
    public void testEncodeDecodeToLength() {
        final AlignmentCollectionHandler handler = new AlignmentCollectionHandler();
        final int[][] fromLists = {{}, {1}, {1, 2, 3, 5, 10}, {1, 2, 3}};
        final int[][] to__Lists = {{}, {1}, {1, 2, 3, 5, 10}, {0, 1, 5}};
        final int numLists = fromLists.length;
        for (int i = 0; i < numLists; i++) {
            final IntList fromLengthList = IntArrayList.wrap(fromLists[i]);
            final IntList toLengthList = IntArrayList.wrap(to__Lists[i]);
            handler.encodeToLength(fromLengthList, toLengthList);
            handler.decodeToLength(fromLengthList, toLengthList);
            assertEquals(IntArrayList.wrap(to__Lists[i]), toLengthList);

        }

    }

    @Test
    public void nullCollection() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        assertEquals(null, codec.encode(null));
    }

    @Test
    public void smallCollection() throws IOException {
        final HybridChunkCodec2 codec = new HybridChunkCodec2();
        final AlignmentCollectionHandler handler = new AlignmentCollectionHandler();
        handler.setDebugLevel(1);
        codec.setHandler(handler);
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
        Alignments.AlignmentCollection.Builder collection = buildCollection(examples, false);


        testRoundTripWithBuiltEntries(codec, builtEntries, collection, false, null);
    }

    @Test
    // enable to test soft clips compression
    public void testRoundTripWithSoftCLips() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = buildCollection(examples, false);

        testRoundTripWithBuiltEntries(codec, builtEntries, collection, false, exampleClips);
    }

    private class SoftClip {
        String left;
        String right;
        private ByteString leftQual;
        private ByteString rightQual;

        private SoftClip(String left, String right, String leftQual, String rightQual) {
            this.left = left;
            this.right = right;
            try {
                this.leftQual = leftQual != null ? ByteString.copyFrom(leftQual, "US-ASCII") : null;
                this.rightQual = rightQual != null ? ByteString.copyFrom(rightQual, "US-ASCII") : null;
            } catch (UnsupportedEncodingException e) {
                fail("Cannot convert strings to quality scores");
            }
        }
    }

    SoftClip[] exampleClips = {
            new SoftClip("AC", null, "aa", null),
            new SoftClip("AACC", null, "abcd", null),
            new SoftClip("AACC", "TCGGGGG", "haha", "#######"),
            new SoftClip("AACC", null, "hhhh", null),
            new SoftClip("AACC", "TCGGGGG", "haha", "#######"),
    };

    private void addSoftClips(SoftClip[] exampleClips, Alignments.AlignmentCollection.Builder collection) {
        int softClipIndex = 0;
        for (int i = 0; i < collection.getAlignmentEntriesCount(); i++) {
            Alignments.AlignmentEntry.Builder element = collection.getAlignmentEntriesBuilder(i);
            if (exampleClips[softClipIndex].left != null) {
                element.setSoftClippedBasesLeft(exampleClips[softClipIndex].left);
            }
            if (exampleClips[softClipIndex].leftQual != null) {
                element.setSoftClippedQualityLeft(exampleClips[softClipIndex].leftQual);
            }
            if (exampleClips[softClipIndex].right != null) {
                element.setSoftClippedBasesRight(exampleClips[softClipIndex].right);
            }
            if (exampleClips[softClipIndex].rightQual != null) {
                element.setSoftClippedQualityRight(exampleClips[softClipIndex].rightQual);
            }
            softClipIndex++;
            if (softClipIndex > exampleClips.length) {
                softClipIndex = 0;
            }
        }
    }

    @Test

    public void testRoundTrip2() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = buildCollection(examples, true);

        testRoundTripWithBuiltEntries(codec, builtEntries, collection, false, null);
    }

    AlignmentExample[] examplesWithDuplicates = {
            new AlignmentExample(0, 1, 33, 1, "TA", "GG", 27, 31, "--", "TA", 3, 10, 0),    // R1
            new AlignmentExample(0, 2, 33, 2, "TA", "GG", 27, 31, "--", "TA", 3, 10, 0),    // R1
            new AlignmentExample(0, 3, 33, 3, "TA", "GG", 27, 31, "--", "TA", 3, 10, 0),    // R1
            new AlignmentExample(0, 15, 33, 4, "TA", "GG", 27, 31, "--", "TA", 3, 10, 0),   // R1
            new AlignmentExample(0, 101, 9, 5, "..", "AA", 42, 24, "T-", "G.", 5, 11, 0),     // R2
            new AlignmentExample(0, 112, 9, 6, "..", "AA", 42, 24, "T-", "G.", 5, 11, 0),    // R2
            new AlignmentExample(0, 112, 9, 7, "..", "AA", 42, 24, "T-", "G.", 5, 11, 0),    // R2
            new AlignmentExample(1, 1111, 9, 8, "AA--AA", "CCCC", 13, 24, "T-", "G.", 5, 3, 0)   //R3
    };

    @Test
    public void testRoundTripMultipleEntriesDuplicate() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = buildCollection(examplesWithDuplicates, false);
        //  System.out.println(collection.build().toString());
        assertRoundTripMatchExpected(codec, collection);
    }

    @Test
    public void roundTripMore() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        final AlignmentCollectionHandler handler = new AlignmentCollectionHandler();
        handler.setDebugLevel(1);
        codec.setHandler(handler);
        Alignments.AlignmentCollection.Builder collection = loadCollection("test-data/seq-var-test/kevin-synth/sorted-seq-var-reads-gsnap.entries");

        assertRoundTripMatchExpected(codec, collection);
    }

    @Test
    public void testCodeDecode() {

        int q = 3;
        int t = 5;
        int code = AlignmentCollectionHandler.modelQueryAlignedLength(q, t);
        int qdecode = AlignmentCollectionHandler.decodeQueryAlignedLength(code, t);
        assertEquals(qdecode, q);
    }

    @Test
    public void testVarPositionDeltaModTransform() {

        final AlignmentCollectionHandler handler = new AlignmentCollectionHandler();
        final int[][] lists = {{3, 1, 2, 3, 5, 1, 30, 3},
                {3, 1}, {1}, {10, 11}, {}};
        for (int[] list : lists) {

            IntArrayList inputList = IntArrayList.wrap(list);
            //System.out.println("Testing "+inputList);
            IntArrayList output = handler.deltaModTransform(inputList);
            IntList backTransformed = new IntArrayList();
            handler.decodeDeltaModTransform(output, backTransformed);
            assertEquals("the decoded list must match the input list for " + inputList, inputList, backTransformed);
        }
    }

      @Test
    public void testVarPositionDeltaModTransform2() {

        final AlignmentCollectionHandler handler = new AlignmentCollectionHandler();
        final int[][] lists = {{49, 59, 31, 49, 42, 10, 18, 37, 0, 44, 54, 39, 43, 20, 24, 9, 29, 5, 13, 17, 33, 54, 36, 31, 11, 30, 4, 14, 5, 54, 60, 38, 1, 4, 57, 18, 74, 31, 66, 73, 17, 31, 57, 63, 9, 54, 57, 60, 26, 42, 5, 76, 47, 10, 24, 10, 31, 34, 5, 26, 49, 24, 20, 26, 34, 21, 9, 21, 3, 23, 7, 57, 12, 43, 13, 20, 34, 25, 4, 6, 44, 51, 22, 55, 4, 30, 63, 76, 3, 3, 23, 66, 7, 1, 2, 27, 2, 6, 1, 3, 25, 1, 68, 16, 28, 19, 34, 62, 8, 20, 50, 56, 19, 57, 24, 32, 60, 12, 47, 29, 65, 35, 21, 30, 26, 30, 39, 17, 30, 26, 42, 65, 18, 51, 8, 15, 19, 55, 68, 1, 3, 21, 19, 41, 15, 63, 24, 6, 28, 64, 49, 38, 13, 35, 40, 9, 55, 31, 3, 1, 42, 46, 58, 12, 17, 22, 52, 39, 53, 1, 6, 29, 13, 59, 26, 48, 1, 13, 3, 14, 20, 58, 16, 34, 45, 51, 13, 10, 18, 23, 46, 62, 1, 7, 27, 54, 15, 28, 8, 24, 32, 11, 44, 50, 20, 29, 55, 52, 28, 16, 67, 65, 64, 64, 63, 62, 62, 61, 61, 61, 61, 60, 46, 58, 57, 9, 66, 56, 14, 52, 52, 33, 48, 46, 70, 45, 44, 4, 41, 40, 9, 25, 1, 38, 36, 30, 34, 33, 32, 68, 71, 74, 68, 66, 66, 29, 26, 28, 11, 20, 28, 28, 25, 49, 2, 12, 37, 55, 8, 61, 23, 60, 51, 21, 57, 19, 54, 13, 12, 34, 12, 49, 33, 3, 11, 45, 5, 13, 43, 18, 11, 35, 35, 22, 26, 34, 33, 32, 31, 32, 34, 25, 23, 22, 28, 23, 32, 33, 35, 64, 1, 5, 30, 62, 1, 1, 31, 27, 3, 30, 2, 13, 11, 5, 15, 11, 33, 30, 4, 4, 14, 33, 3, 21, 29, 23, 70, 67, 6, 2, 1, 5, 9, 47, 67, 13, 9, 43, 49, 1, 41, 38, 37, 35, 21, 6, 30, 29, 32, 24, 18, 23, 8, 7, 9, 18, 23, 4, 25, 2, 1, 22, 67, 20, 2, 20, 23, 43, 58, 60, 31, 53, 73, 69, 48, 19, 68, 23, 63, 66, 44, 6, 60, 61, 58, 8, 15, 16, 33, 42, 56, 14, 53, 50, 47, 59, 66, 5, 17, 46, 23, 29, 41, 48, 19, 31, 29, 29, 28, 42, 27, 27, 43, 25, 11, 19, 34, 3, 16, 21, 21, 21, 19, 24, 28, 18, 17, 17, 20, 16, 12, 36, 9, 10, 6, 3, 10, 12, 3, 45, 8, 5, 19, 1, 4, 29, 12, 7, 16, 38, 8, 14, 3, 68, 3, 6, 24, 1, 65, 65, 25, 3, 1, 38, 33, 66, 1, 15, 17, 11, 33, 57, 41, 15, 20, 13, 20, 41, 25, 63, 1, 9, 7, 51, 19, 14, 45, 24, 20, 45, 1, 61, 1, 12, 58, 61, 15, 15, 31, 41, 47, 29, 18, 30, 61, 41, 9, 8, 12, 24, 4, 7, 13, 53, 7, 27, 9, 16, 27, 1, 5, 23, 17, 13, 20, 25, 26, 40, 14, 21, 13, 14, 20, 27, 13, 38, 37, 33, 11, 30, 3, 11, 8, 34, 1, 4, 1, 59, 64, 66, 24, 31, 20, 44, 3, 18, 59, 22, 14, 27, 34, 38, 9, 29, 72, 75, 47, 38, 43, 45, 61, 30, 15, 15, 30, 2, 15, 3, 25, 4, 33, 63, 35, 27, 23, 18, 4, 37, 34, 60, 25, 15, 8, 1, 41, 1, 6, 28, 16, 50, 60, 12, 23, 33, 8, 12, 57, 7, 32, 23, 32, 42, 62, 17, 20, 28, 10, 76, 1, 49, 1, 6, 33, 54, 50, 17, 23, 26, 2, 47, 3, 18, 46, 4, 29, 46, 67, 3, 45, 63, 46, 63, 66, 24, 61, 25, 57, 31, 47, 20, 26, 14, 33, 50, 37, 66, 48, 74, 44, 27, 4, 11, 1, 9, 73, 75, 29, 30, 34, 1, 1, 67, 13, 37, 24, 6, 8, 10, 1, 1, 25, 55, 63, 33, 3, 17, 21, 19, 60, 27, 17, 55, 53, 66, 32, 33, 38, 39, 23, 47, 57, 37, 15, 55, 68, 26, 14, 39, 31, 53, 32, 42, 5, 9, 68, 5, 28, 4, 1, 19, 3, 27, 40, 6, 34, 13, 47, 5, 14, 47, 54, 59, 65, 23, 30, 57, 34, 19, 25, 17, 6, 63, 10, 2, 5, 47, 41, 1, 67, 51, 18, 59, 19, 24, 29, 35, 61, 45, 61, 50, 56, 60, 16, 20, 28, 57, 34, 35, 18, 17, 55, 68, 72, 74, 4, 47, 17, 17, 68, 59, 10, 46, 8, 33, 47, 32, 7, 52, 43, 31, 14, 16, 1, 4, 45, 64, 3, 13, 1, 2, 4, 1, 6, 32, 49, 22, 57, 18, 16, 31, 32, 67, 25, 65, 76, 1, 7, 15, 5, 34, 30, 26, 18, 28, 20, 28, 34, 33, 40, 3, 34, 21, 54, 66, 48, 33, 7, 7, 2, 64, 55, 23, 17, 16, 61, 54, 61, 49, 55, 7, 9, 14, 44, 23, 33, 19, 10, 10, 51, 22, 6, 64, 61, 10, 22, 26, 33, 48, 15, 36, 2, 45, 32, 58, 25, 31, 26, 21, 4, 18, 40, 3, 67, 63, 17, 19, 66, 21, 1, 39, 66, 4, 58, 66, 51, 53, 66, 55, 48, 11, 17, 49, 19, 63, 39, 9, 41, 51, 56, 3, 10, 30, 50, 25, 7, 4, 1, 33, 1, 8, 6, 12, 41, 67, 1, 27, 29, 21, 1, 23, 27, 38, 9, 12, 14, 4, 64, 66, 53, 67, 73, 11, 19, 7, 28, 68, 24, 2, 4, 68, 71, 74, 7, 68, 3, 9, 1, 35, 47, 26, 25, 3, 11, 12, 30, 33, 14, 58, 51, 46, 6, 12, 51, 11, 19, 32, 6, 56, 4, 15, 26, 35, 37, 50, 25, 51, 16, 22, 42, 7, 1, 6, 1, 39, 40, 5, 44, 23, 48, 48, 31, 75, 28, 32, 41, 17, 26, 50, 13, 25, 14, 8, 32, 49, 5, 52, 68, 3, 58, 11, 15, 21, 27, 49, 30, 14, 16, 68, 69, 26, 1, 29, 1, 20, 20, 26, 3, 40, 2, 17, 49, 17, 3, 1, 27, 18, 20, 48, 3, 31, 60, 49, 14, 15, 2, 9, 1, 1, 16, 28, 12, 17, 23, 3, 27, 13, 45, 48, 13, 1, 26, 42, 26, 23, 29, 31, 8, 32, 23, 56, 19, 6, 12, 29, 33, 32, 2, 30, 27, 25, 17, 1, 2, 28, 31, 57, 5, 41, 16, 4, 32, 26, 5, 32, 14, 12, 56, 52, 57, 37, 13, 47, 34, 15, 48, 30, 34, 4, 74, 52, 1, 43, 5, 6, 17, 1, 9, 60, 16, 34, 54, 27, 60, 17, 55, 22, 65, 34, 41, 46, 4, 20, 51, 8, 51, 12, 6, 67, 5, 7, 33, 42, 1, 39, 47, 59, 63, 1, 40, 4, 28, 30, 17, 59, 1, 27, 7, 19, 36, 57, 25, 57, 28, 43, 15, 3, 43, 4, 18, 16, 10, 2, 15, 72, 32, 15, 14, 34, 4, 43, 67, 33, 57, 63, 14, 16, 45, 30, 9, 36, 46, 54, 58, 6, 35, 65, 49, 34, 17, 9, 20, 28, 22, 11, 1, 47, 68, 14, 20, 1, 10, 10, 1, 30, 53, 5, 1, 58, 28, 15, 28, 33, 29, 5, 35, 40, 6, 53, 14, 34, 2, 18, 33, 4, 21, 62, 25, 29, 37, 34, 16, 15, 24, 22, 47, 54, 18, 5, 17, 63, 66, 13, 29, 67, 23, 25, 18, 9, 33, 56, 32, 13, 3, 26, 55, 48, 12, 26, 13, 3, 40, 15, 56, 7, 26, 61, 20, 34, 11, 37, 54, 13, 29, 26, 42, 56, 1, 3, 24, 14, 61, 26, 28, 33, 20, 24, 13, 20, 34, 45, 10, 65, 9, 64, 1, 61, 14, 40, 25, 30, 21, 24, 18, 65, 20, 10, 7, 4, 43, 13, 24, 23, 58, 4, 65, 11, 11, 57, 64, 61, 30, 7, 58, 31, 44, 50, 61, 12, 49, 35, 50, 6, 54, 66, 1, 50, 62, 1, 9, 33, 36, 9, 30, 42, 6, 31, 42, 13, 50, 30, 6, 65, 22, 32, 60, 47, 33, 17, 41, 37, 65, 11, 47, 21, 25, 34, 17, 31, 15, 4, 34, 65, 12, 11, 61, 76, 47, 3, 8, 15, 21, 28, 68, 64, 13, 45, 12, 15, 26, 5, 13, 16, 25, 57, 31, 48, 41, 5, 23, 28, 36, 13, 18, 26, 26, 7, 35, 5, 8, 18, 26, 28, 32, 29, 46, 21, 40, 63, 49, 61, 49, 58, 48, 60, 12, 28, 33, 22, 34, 58, 64, 13, 25, 36, 4, 6, 13, 58, 15, 35, 3, 7, 17, 40, 33, 27, 14, 56, 58, 12, 54, 56, 37, 39, 21, 23, 31, 46, 57, 12, 13, 12, 13, 49, 34, 23, 26, 11, 44, 27, 34, 25, 14, 37, 47, 14, 17, 24, 29, 57, 8, 18, 23, 31, 15, 17, 29, 6, 19, 66, 17, 20, 52, 61, 49, 9, 24, 38, 41, 51, 6, 21, 35, 48, 15, 19, 23, 32, 13, 24, 12, 23, 5, 43, 54, 35, 46, 28, 39, 58, 7, 37, 56, 29, 32, 35, 1, 9, 28, 71, 16, 33, 13, 65, 73, 49, 46, 54, 39, 36, 44, 35, 43, 34, 22, 30, 66, 3, 5, 2, 36, 10, 59, 4, 3, 52, 38, 36, 3, 21, 13, 31, 36, 47, 51, 52, 58, 54, 39, 51, 22, 49, 46, 11, 23, 61, 22, 59, 21, 59, 16, 44, 25, 54, 8, 5, 16, 18, 8, 17, 25, 36, 15, 17, 34, 31, 38, 3, 12, 2, 18, 7, 36, 23, 58, 1, 3, 1, 1, 48, 1, 48, 1, 47, 8, 14, 19, 22, 49, 49, 7, 22, 25, 15, 3, 2, 15, 3, 20, 35, 3, 46, 36, 5, 11, 17, 51, 11, 1, 3, 42, 49, 3, 22, 24, 55, 62, 21, 29, 30, 5, 68, 31, 10, 56, 49, 2, 27, 10, 13, 18, 22, 28, 16, 20, 22, 38, 2, 18, 25, 20, 27, 20, 27, 33, 2, 9, 1, 38, 11, 1, 67, 28, 20, 37, 1, 8, 36, 28, 55, 1, 33, 28, 49, 52, 68, 66, 16, 23, 21, 26, 60, 59, 67, 31, 55, 58, 56, 4, 12, 17, 21, 3, 16, 29, 19, 12, 16, 66, 2, 5, 25, 35, 54, 22, 70, 68, 19, 67, 59, 1, 59, 55, 54, 53, 48, 44, 43, 40, 40, 25, 25, 24, 16, 18, 30, 15, 13, 24, 1, 17, 26, 11, 35, 9, 16, 61, 7, 26, 60, 4, 62, 55, 20, 24, 33, 41, 22, 50, 53, 66, 14, 38, 13, 37, 11, 35, 6, 10, 15, 18, 23, 19, 2, 2, 16, 14, 13, 13, 12, 12, 11, 9, 9, 2, 8, 5, 8, 68, 37, 13, 1, 7, 68, 1, 67, 28, 62, 61, 59, 59, 61, 56, 58, 61, 11, 13, 51, 54, 56, 54, 56, 51, 47, 49, 42, 41, 43, 38, 40, 13, 31, 33, 27, 29, 8, 12, 25, 24, 26, 30, 21, 21, 13, 16, 18, 20, 26, 14, 16, 18, 24, 13, 15, 12, 58, 10, 9, 11, 1, 7, 9, 7, 9, 5, 7, 5, 7, 3, 1, 3, 37, 4, 7, 55, 39, 55, 35, 51, 50, 14, 48, 33, 36, 37, 12, 16, 9, 6, 6, 24, 27, 23, 5, 16, 19, 1, 15, 1, 20, 68, 1, 10, 15, 65, 6, 23, 25, 29, 4, 53, 58, 75, 70, 53, 52, 18, 48, 26, 32, 22, 24, 13, 50, 14, 14, 19, 5, 3, 4, 8, 8, 5, 16, 33, 1, 15, 42, 32, 42, 46, 68, 31, 48, 66, 65, 64, 44, 56, 54, 54, 53, 62, 13, 48, 43, 38, 37, 36, 31, 29, 41, 27, 24, 24, 46, 9, 19, 22, 2, 26, 17, 3, 34, 52, 2, 12, 1, 14, 23, 11, 20, 5, 12, 20, 6, 16, 19, 5, 7, 10, 18, 5, 17, 3, 10, 16, 3, 16, 4, 7, 15, 8, 10, 13, 1, 12, 10, 10, 8, 16, 18, 3, 6, 6, 2, 1, 1, 1, 14, 7, 12, 6, 3, 1, 23, 3, 9, 3, 14, 3, 3, 14, 33, 1, 50, 15, 1, 3, 20, 47, 3, 14, 3, 14, 45, 21, 41, 17, 44, 2, 44, 7, 4, 7, 11, 23, 3, 10, 16, 20, 24, 44, 51, 66, 40, 44, 30, 5, 43, 39, 47, 34, 42, 25, 31, 29, 37, 23, 31, 22, 30, 3, 17, 20, 57, 38, 50, 18, 25, 46, 29, 41, 41, 62, 37, 52, 65, 23, 35, 22, 37, 41, 22, 37, 50, 21, 36, 49, 12, 20, 27, 38, 51, 10, 25, 38, 11, 7, 18, 31, 13, 21, 13, 17, 38, 1, 55, 8, 16, 4, 12, 50, 1, 25, 47, 56, 11, 44, 10, 15, 19, 27, 37, 62, 3, 11, 21, 47, 48, 65, 32, 42, 26, 30, 21, 14, 7, 6, 16, 19, 3, 46, 58, 49, 34, 20, 22, 26, 6, 3, 32, 13, 6, 4, 11, 23, 18, 24, 29, 64, 24, 1, 56, 1, 42, 5, 8, 14, 38, 15, 4, 30, 5, 31, 40, 43, 49, 54, 51, 4, 8, 34, 45, 49, 52, 57, 27, 41, 13, 29, 30, 36, 38, 47, 49, 52, 56, 58, 49, 7, 5, 10, 45, 1, 11, 16, 39, 33, 16, 26, 40, 51, 4, 17, 24, 29, 53, 9, 11, 50, 53, 49, 64, 5, 24, 26, 44, 19, 21, 39, 16, 35, 50, 5, 8, 11, 13, 1, 57, 43, 8, 13, 16, 18, 20, 22, 25, 47, 29, 33, 41, 31, 16, 28, 32, 13, 10, 48, 30, 34, 47, 35, 25, 35, 27, 46, 4, 9, 14, 16, 23, 22, 28, 31, 1, 24, 32, 31, 13, 17, 22, 33, 57, 29, 35, 34, 6, 17, 9, 22, 24, 30, 23, 28, 19, 31, 42, 42, 5, 17, 63, 21, 20, 26, 26, 34, 36, 45, 1, 13, 63, 28, 34, 18, 26, 68, 17, 23, 27, 1, 10, 8, 12, 18, 13, 31, 39, 20, 50, 56, 71, 24, 33, 4, 44, 3, 59, 25, 27, 58, 40, 55, 15, 4, 4, 1, 3, 3, 2, 2, 1, 1, 1, 1, 62, 11, 17, 55, 59, 37, 15, 28, 58, 43, 19, 53, 10, 31, 11, 15, 19, 28, 2, 5, 72, 75, 16, 9, 11, 16, 25, 52, 1, 22, 18, 35, 7, 9, 66, 12, 20, 31, 34, 34, 9, 62, 15, 60, 51, 9, 17, 27, 21, 43, 68, 19, 39, 49, 55, 67, 31, 20, 30, 2, 8, 6, 8, 17, 28, 4, 9, 15, 18, 24, 27, 9, 5, 12, 15, 19, 6, 7, 57, 4, 25, 1, 30, 7, 18, 4, 9, 30, 9, 55, 69, 8, 9, 15, 22, 25, 28, 15, 19, 25, 12, 9, 36, 4, 49, 31, 16, 12, 52, 1, 68, 62, 41, 56, 18, 32, 39, 1, 13, 25, 28, 24, 22, 55, 9, 4, 4, 20, 8, 23, 7, 63, 19, 12, 39, 9, 41, 63, 3, 11, 63, 68, 64, 9, 11, 24, 20, 57, 22, 1, 28, 66, 1, 44, 50, 27, 8, 1, 30, 5, 13, 38, 3, 5, 23, 35, 1, 9, 6, 32, 20, 2, 15, 54, 4, 28, 31, 26, 25, 30, 18, 5, 10, 14, 16, 24, 47, 67, 11, 16, 31, 8, 48, 30, 22, 62, 34, 5, 29, 51, 57, 46, 3, 21, 26, 51, 16, 22, 27, 8, 10, 51, 68, 19, 29, 20, 4, 13, 50, 14, 61, 33, 1, 20, 22, 6, 61, 67, 26, 25, 40, 32, 40, 59, 5, 11, 15, 18, 20, 16, 68, 66, 23, 5, 18, 22, 53, 3, 24, 31, 18, 30, 28, 25, 13, 25, 24, 5, 21, 29, 12, 11, 11, 1, 1, 1, 25, 3, 17, 31, 17, 24, 37, 1, 27, 28, 30, 43, 6, 18, 21, 56, 5, 11, 50, 12, 49, 64, 1, 17, 1, 30, 3, 8, 12, 5, 33, 8, 3, 15, 33, 4, 3, 25, 28, 30, 7, 17, 27, 4, 8, 14, 20, 24, 41, 54, 63, 67, 9, 19, 8, 12, 5, 34, 3, 8, 7, 35, 8, 12, 23, 23, 27, 8, 14, 13, 28, 31, 28, 31, 2, 22, 12, 1, 31, 50, 76, 24, 25, 1, 68, 20, 66, 9, 6, 35, 15, 26, 15, 31, 6, 75, 74, 36, 49, 68, 31, 54, 33, 52, 16, 26, 12, 17, 5, 15, 19, 22, 32, 4, 11, 16, 18, 44, 43, 43, 13, 26, 34, 24, 32, 31, 31, 30, 27, 34, 32, 33, 26, 31, 23, 31, 27, 22, 27, 5, 8, 21, 21, 25, 19, 25, 31, 20, 35, 17, 19, 10, 9, 16, 18, 20, 29, 31, 10, 20, 26, 29, 10, 14, 16, 18, 21, 25, 29, 31, 44, 34, 15, 14, 27, 5, 15, 21, 53, 68, 19, 50, 56, 63, 65, 44, 33, 14, 21, 32, 18, 9, 12, 49, 19, 9, 6, 12, 47, 7, 37, 4, 23, 60, 67, 34, 25, 76, 40, 68, 14, 31, 9, 11, 62, 23, 30, 50, 30, 33, 4, 9, 14, 16, 23, 4, 25, 60, 59, 37, 1, 17, 28, 22, 39, 16, 13, 19, 17, 29, 1, 50, 14, 6, 23, 13, 21, 48, 32, 54, 15, 18, 7, 30, 1, 5, 7, 11, 20, 37, 1, 18, 26, 30, 50, 8, 14, 12, 14, 3, 17, 1, 1, 35, 34, 35, 35, 9, 29, 7, 22, 27, 32, 17, 26, 26, 65, 10, 19, 31, 5, 7, 14, 26, 8, 11, 32, 34, 40, 21, 16, 48, 54, 39, 4, 16, 11, 16, 6, 2, 3, 22, 7, 27, 11, 14, 44, 6, 13, 1, 37, 23, 3, 12, 18, 20, 4, 3, 1, 3, 4, 8, 14, 18, 38, 15, 26, 50, 18, 61, 1, 15, 43, 54, 17, 49, 55, 31, 5, 34, 19, 19, 33, 32, 35, 35, 14, 31, 28, 8, 11, 26, 67, 19, 21, 14, 11, 16, 27, 34, 6, 8, 21, 11, 23, 28, 3, 26, 23, 68, 9, 11, 22, 46, 9, 27, 44, 32, 35, 59, 58, 1, 3, 15, 4, 6, 10, 16, 19, 23, 15, 6, 62, 1, 34, 8, 24, 12, 27, 2, 15, 6, 10, 13, 17, 23, 14, 26, 30, 33, 24, 59, 65, 64, 31, 48, 3, 19, 38, 42, 3, 11, 51, 51, 27, 24, 13, 19, 23, 28, 32, 15, 16, 19, 24, 39, 17, 21, 23, 53, 27, 6, 68, 67, 66, 64, 3, 63, 61, 60, 60, 29, 54, 57, 49, 64, 55, 50, 49, 49, 52, 48, 47, 46, 44, 44, 43, 42, 40, 18, 15, 28, 28, 28, 25, 3, 20, 31, 16, 14, 10, 12, 1, 20, 6, 10, 18, 20, 31, 68, 55, 44, 41, 32, 32, 24, 24, 14, 8, 10, 37, 31, 15, 17, 28, 1, 12, 28, 21, 16, 49, 62, 4, 32, 29, 1, 6, 68, 54, 5, 10, 39, 11, 32, 5, 9, 26, 35, 31, 36, 57, 8, 32, 56, 10, 22, 37, 1, 14, 18, 6, 23, 4, 7, 37, 41, 41, 59, 41, 68, 48, 60, 66, 10, 34, 22, 25, 24, 10, 23, 29, 3, 14, 1, 24, 10, 17, 20, 23, 25, 30, 35, 38, 9, 35, 37, 16, 9, 10, 12, 31, 1, 66, 13, 24, 47, 23, 51, 10, 15, 20, 23, 25, 29, 12, 5, 18, 46, 6, 18, 26, 30, 51, 68, 58, 5, 8, 13, 17, 47, 25, 9, 3, 24, 17, 35, 14, 6, 35, 40, 56, 58, 65, 65, 41, 28, 33, 18, 19, 23, 30, 4, 8, 18, 30, 11, 28, 26, 32, 22, 27, 33, 7, 21, 33, 4, 33, 32, 7, 10, 12, 28, 10, 13, 21, 27, 29, 42, 64, 13, 29, 22, 23, 10, 22, 30, 9, 17, 23, 28, 3, 14, 27, 10, 18, 30, 4, 7, 10, 17, 22, 13, 3, 16, 24, 13, 23, 1, 5, 15, 17, 31, 19, 23, 31, 38, 46, 35, 22, 29, 32, 2, 8, 2, 67, 10, 18, 49, 31, 7, 4, 31, 8, 7, 12, 21, 11, 20, 23, 30, 32, 37, 34, 6, 14, 20, 26, 28, 52, 1, 3, 31, 5, 49, 69, 5, 57, 19, 11, 15, 19, 31, 45, 1, 5, 31, 27, 47, 49, 41, 21, 25, 32, 11, 16, 21, 11, 16, 18, 21, 23, 30, 33, 58, 60, 36, 36, 22, 30, 68, 7, 12, 16, 20, 27, 29, 11, 20, 25, 22, 7, 61, 22, 16, 21, 26, 44, 23, 50, 7, 60, 38, 40, 8, 42, 28, 39, 2, 10, 1, 4, 16, 25, 40, 45, 48, 8, 11, 13, 16, 22, 31, 23, 2, 32, 2, 26, 7, 6, 21, 65, 23, 7, 64, 35, 31, 34, 44, 1, 1, 37, 68, 71, 37, 14, 29, 1, 74, 76, 10, 53, 23, 14, 49, 5, 13, 18, 20, 19, 21, 6, 11, 14, 40, 21, 64, 31, 41, 39, 56, 9, 17, 19, 23, 3, 6, 11, 21, 24, 6, 1, 48, 47, 14, 31, 12, 30, 28, 11, 1, 6, 3, 3, 2, 68, 24, 13, 40, 4, 48, 49, 63, 58, 33, 27, 33, 19, 34, 6, 10, 16, 22, 28, 6, 11, 1, 37, 4, 59, 28, 29, 33, 25, 4, 68, 12, 60, 49, 43, 42, 15, 17, 26, 13, 21, 58, 55, 31, 12, 2, 1, 67, 26, 18, 22, 24, 27, 10, 62, 27, 5, 22, 52, 29, 55, 25, 6, 5, 21, 11, 5, 50, 55, 30, 16, 15, 5, 21, 35, 59, 66, 47, 48, 50, 66, 1, 6, 1, 11, 9, 67, 21, 3, 11, 14, 4, 33, 1, 8, 1, 6, 50, 4, 1, 28, 32, 33, 43, 60, 4, 40, 12, 20, 3, 8, 50, 14, 17, 13, 29, 33, 46, 39, 35, 28, 27, 24, 30, 29, 7, 11, 17, 27, 23, 22, 66, 21, 21, 21, 17, 13, 12, 12, 9, 7, 7, 5, 2, 1, 58, 3, 25, 1, 8, 6, 29, 17, 1, 19, 29, 35, 27, 58, 12, 35, 3, 15, 22, 26, 29, 48, 9, 6, 65, 72, 1, 2, 8, 40, 11, 58, 5, 61, 55, 31, 41, 21, 3, 9, 17, 34, 45, 14, 14, 34, 66, 59, 18, 6, 11, 24, 1, 13, 15, 43, 34, 29, 3, 34, 35, 31, 7, 22, 30, 33, 25, 35, 10, 22, 15, 40, 11, 38, 3, 18, 21, 35, 42, 58, 4, 7, 13, 22, 24, 27, 22, 16, 21, 47, 57, 16, 20, 25, 27, 30, 13, 29, 39, 4, 6, 19, 57, 54, 37, 31, 25, 1, 30, 13, 25, 3, 11, 14, 22, 30, 16, 40, 53, 19, 51, 1, 14, 26, 63, 16, 12, 10, 65, 6, 37, 34, 15, 14, 17, 64, 27, 11, 2, 60, 55, 42, 61, 34, 61, 68, 8, 17, 16, 23, 1, 2, 50, 47, 35, 24, 3, 10, 14, 5, 11, 5, 7, 11, 13, 15, 6, 13, 15, 8, 36, 41, 26, 14, 34, 4, 66, 43, 28, 12, 7, 20, 50, 19, 55, 37, 41, 11, 22, 1, 3, 5, 10, 20, 22, 24, 9, 30, 30, 28, 20, 43, 60, 68, 1, 13, 43, 32, 8, 13, 26, 27, 42, 37, 30, 1, 42, 1, 1, 4, 10, 5, 10, 1, 14, 21, 29, 31, 44, 5, 8, 19, 22, 29, 3, 13, 49, 11, 17, 27, 20, 34, 7, 16, 20, 31, 16, 3, 8, 17, 32, 7, 19, 40, 29, 31, 28, 33, 70, 47, 34, 41, 45, 49, 48, 12, 27, 68, 21, 11, 35, 16, 31, 55, 9, 27, 33, 1, 31, 33, 57, 5, 36, 23, 36, 4, 11, 19, 25, 20, 14, 12, 26, 32, 67, 9, 36, 3, 11, 44, 64, 68, 3, 8, 31, 49, 15, 26, 39, 46, 50, 10, 19, 23, 4, 11, 22, 26, 1, 1, 11, 11, 35, 4, 34, 61, 69, 5, 13, 26, 3, 5, 20, 40, 39, 45, 59, 18, 22, 34, 14, 1, 63, 44, 31, 71, 24, 61, 14, 61, 58, 36, 34, 69, 10, 15, 39, 65, 25, 29, 28, 46, 19, 22, 24, 47, 41, 37, 4, 12, 52, 16, 47, 39, 48, 4, 13, 1, 10, 29, 37, 7, 31, 7, 11, 14, 16, 24, 47, 45, 44, 24, 61, 8, 20, 27, 31, 29, 4, 62, 21, 11, 19, 33, 19, 35, 66, 55, 6, 14, 65, 63, 59, 57, 50, 57, 56, 50, 50, 48, 52, 27, 17, 26, 21, 21, 20, 16, 16, 15, 12, 5, 5, 28, 39, 13, 33, 12, 30, 28, 28, 7, 28, 1, 12, 14, 17, 31, 24, 29, 26, 31, 24, 36, 21, 18, 11, 9, 29, 40, 62, 42, 42, 50, 7, 20, 31, 40, 53, 66, 3, 35, 40, 58, 2, 42, 45, 1, 74, 30, 35, 31, 21, 26, 45, 22, 22, 45, 22, 53, 67, 13, 41, 49, 44, 45, 22, 4, 6, 14, 20, 30, 4, 6, 7, 12, 63, 66, 30, 1, 19, 14, 44, 12, 66, 62, 56, 47, 13, 16, 34, 52, 40, 14, 69, 39, 9, 56, 43, 48, 3, 32, 32, 48, 17, 47, 18, 15, 12, 39, 3, 2, 66, 1, 60, 57, 56, 51, 51, 49, 49, 44, 41, 67, 41, 67, 40, 66, 33, 59, 32, 58, 30, 56, 29, 55, 29, 55, 23, 49, 20, 46, 20, 46, 14, 40, 5, 31, 15, 7, 7, 55, 7, 6, 4, 50, 3, 31, 14, 29, 31, 26, 23, 7, 26, 5, 30, 35, 19, 60, 1, 5, 13, 27, 49, 34, 33, 58, 3, 25, 25, 20, 17, 15, 15, 13, 13, 9, 56, 5, 45, 68, 34, 68, 50, 68, 68, 63, 19, 62, 58, 36, 34, 25, 14, 23, 5, 17, 5, 16, 9, 8, 2, 42, 1, 41, 40, 38, 38, 37, 35, 27, 27, 22, 18, 12, 10, 37, 30, 5, 51, 21, 25, 63, 53, 44, 49, 43, 35, 24, 23, 17, 13, 10, 6, 1, 9, 23, 62, 35, 70, 72, 68, 29, 67, 65, 67, 74, 63, 65, 23, 34, 49, 51, 58, 45, 47, 54, 44, 46, 53, 67, 31, 33, 40, 54, 66, 28, 30, 37, 51, 28, 30, 37, 51, 63, 26, 28, 35, 23, 25, 32, 46, 58, 21, 23, 30, 44, 56, 17, 19, 23, 26, 15, 17, 24, 38, 50, 12, 14, 21, 35, 12, 14, 21, 35, 47, 7, 9, 16, 30, 42, 4, 6, 13, 7, 21, 33, 7, 21, 33, 7, 21, 33, 59, 7, 21, 33, 4, 18, 30, 56, 60, 3, 17, 29, 55, 59, 8, 12, 24, 9, 21, 47, 51, 8, 20, 7, 19, 45, 49, 6, 18, 44, 48, 5, 17, 43, 47, 4, 16, 42, 46, 4, 16, 42, 46, 3, 15, 41, 45, 1, 13, 39, 43, 1, 13, 39, 43, 10, 36, 40, 66, 4, 30, 34, 60, 4, 30, 34, 60, 24, 28, 54, 19, 23, 42, 49, 13, 17, 43, 11, 15, 41, 9, 13, 39, 9, 13, 39, 8, 12, 38, 6, 10, 18, 36, 1, 5, 31, 37, 4, 4, 30, 26, 22, 21, 17, 14, 14, 12, 10, 7, 5, 4, 2, 14, 2, 4, 3, 3, 67, 57, 57, 56, 56, 55, 54, 53, 53, 48, 43, 38, 34, 33, 3, 17, 28, 31, 29, 28, 26, 25, 22, 19, 18, 17, 16, 16, 15, 15, 14, 8, 6, 5, 5, 2, 2, 15, 31, 50, 23, 31, 33, 68, 67, 62, 63, 52, 48, 51, 46, 39, 39, 35, 34, 29, 29, 60, 29, 28, 26, 25, 24, 17, 14, 13, 74, 13, 74, 52, 52, 52, 52, 41, 19, 39, 39, 38, 38, 36, 31, 42, 29, 29, 28, 27, 26, 25, 9, 11, 22, 24, 28, 17, 14, 14, 13, 13, 12, 7, 68, 19, 63, 56, 4, 26, 1, 65, 7, 10, 20, 63, 19, 22, 1, 4, 43, 60, 5, 1, 65, 49, 40, 7, 51, 42, 12, 22, 60, 45, 39, 31, 65, 31, 4, 26, 21, 4, 1, 3, 13, 21, 31, 60, 2, 34, 18, 19, 65, 53, 43, 20, 9, 8, 21, 19, 68, 2, 7, 11, 3, 12, 20, 25, 28, 25, 35, 4, 11, 45, 47, 39, 37, 13, 5, 33, 35, 50, 53, 4, 27, 27, 65, 49, 29, 45, 35, 28, 13, 60, 44, 60, 13, 29, 25, 3, 37, 51, 67, 16, 30, 46, 60, 43, 57, 27, 13, 33, 41, 60, 10, 7, 13, 6, 24, 4, 66, 35, 5, 11, 18, 23, 26, 22, 32, 17, 21, 26, 12, 11, 21, 6, 23, 25, 26, 42, 29, 24, 50, 21, 19, 8, 12, 17, 19, 22, 29, 31, 4, 12, 37, 5, 12, 17, 21, 24, 62, 39, 54, 50, 16, 4, 18, 22, 25, 64, 35, 15, 62, 1, 53, 52, 50, 50, 4, 6, 37, 24, 11, 6, 37, 48, 48, 17, 31, 12, 14, 19, 40, 33, 12, 1, 46, 11, 25, 43, 65, 28, 65, 8, 50, 23, 47, 33, 1, 65, 30, 1, 53, 28, 28, 34, 12, 34, 9, 5, 18, 23, 30, 32, 2, 62, 62, 53, 40, 20, 4, 13, 17, 21, 62, 1, 21, 63, 54, 1, 32, 21, 21, 12, 24, 57, 16, 46, 64, 5, 7, 12, 36, 52, 1, 66, 44, 47, 58, 30, 47, 43, 3, 1, 9, 23, 1, 8, 58, 35, 76, 30, 34, 24, 45, 57, 48, 9, 31, 27, 32, 35, 13, 23, 50, 17, 16, 29, 14, 34, 41, 60, 1, 25, 39, 47, 34, 12, 19, 49, 20, 30, 5, 25, 28, 14, 33, 16, 18, 20, 23, 23, 43, 55, 7, 16, 65, 6, 16, 6, 23, 27, 56, 15, 41, 4, 68, 13, 65, 9, 12, 15, 18, 27, 5, 2, 50, 25, 24, 30, 55, 12, 32, 22, 24, 18, 20, 8, 16, 36, 15, 33, 7, 24, 26, 60, 29, 32, 38, 64, 52, 61, 44, 29, 44, 11, 14, 19, 68, 33, 45, 40, 3, 30, 32, 39, 12, 2, 57, 25, 34, 32, 52, 24, 16, 68, 8, 26, 30, 31, 5, 20, 4, 9, 11, 35, 18, 26, 24, 34, 65, 20, 24, 10, 23, 34, 61, 3, 9, 12, 30, 19, 28, 3, 33, 36, 8, 15, 17, 32, 26, 17, 22, 8, 32, 39, 10, 16, 15, 18, 20, 25, 15, 3, 27, 8, 2, 29, 4, 10, 4, 28, 25, 30, 6, 3, 1, 6, 9, 28, 17, 29, 31, 44, 3, 25, 19, 10, 7, 21, 3, 17, 9, 20, 22, 22, 45, 4, 7, 21, 1, 1, 11, 8, 11, 27, 1, 32, 38, 28, 14, 51, 50, 13, 48, 12, 5, 5, 24, 5, 19, 21, 30, 33, 35, 39, 47, 55, 2, 14, 31, 64, 61, 3, 16, 35, 11, 68, 9, 8, 60, 21, 28, 12, 14, 20, 23, 6, 29, 75, 25, 23, 53, 60, 65, 30, 46, 58, 5, 8, 3, 21, 17, 32, 1, 7, 1, 9, 6, 59, 5, 7, 18, 29, 24, 29, 35, 9, 45, 56, 3, 9, 6, 19, 45, 2, 25, 9, 22, 6, 3, 18, 4, 4, 1, 16, 2, 11, 35, 33, 12, 1, 12, 6, 39, 18, 33, 46, 17, 43, 15, 41, 28, 13, 1, 13, 13, 3, 12, 14, 15, 22, 25, 4, 9, 5, 43, 65, 32, 13, 73, 34, 13, 1, 33, 54, 1, 55, 22, 20, 29, 6, 42, 13, 33, 60, 4, 14, 34, 39, 8, 13, 29, 34, 1, 3, 23, 25, 28, 27, 46, 2, 21, 72, 75, 26, 76, 17, 2, 48, 6, 35, 15, 50, 31, 42, 12, 6, 10, 12, 24, 31, 34, 5, 12, 7, 23, 10, 22, 19, 17, 65, 48, 13, 5, 59, 25, 15, 62, 55, 4, 17, 3, 7, 10, 15, 18, 3, 7, 12, 15, 18, 33, 10, 19, 23, 31, 8, 8, 46, 42, 15, 48, 24, 56, 17, 27, 20, 48, 5, 64, 3, 66, 8, 13, 19, 43, 60, 54, 1, 19, 32, 40, 7, 25, 3, 5, 34, 62, 28, 19, 1, 20, 8, 23, 20, 5, 9, 31, 60, 24, 7, 11, 10, 38, 2, 6, 45, 4, 43, 50, 44, 60, 11, 15, 23, 65, 35, 27, 37, 3, 6, 39, 68, 58, 43, 29, 16, 21, 23, 32, 51, 24, 4, 19, 9, 12, 3, 6, 21, 39, 41, 37, 3, 7, 4, 44, 19, 40, 36, 43, 41, 48, 65, 22, 5, 23, 5, 18, 22, 54, 53, 11, 9, 10, 1, 27, 44, 12, 25, 22, 39, 49, 35, 23, 35, 9, 23, 33, 68, 7, 23, 31, 67, 5, 70, 48, 1, 29, 23, 10, 56, 68, 5, 13, 17, 38, 43, 62, 36, 17, 5, 14, 23, 30, 50, 16, 43, 8, 21, 25, 31, 42, 32, 30, 3, 8, 46, 32, 18, 38, 57, 66, 14, 38, 26, 35, 32, 35, 62, 40, 54, 11, 14, 22, 24, 26, 30, 1, 2, 28, 30, 19, 6, 19, 6, 15, 50, 60, 19, 30, 16, 18, 4, 1, 4, 1, 4, 32, 4, 8, 11, 14, 24, 2, 3, 14, 6, 23, 3, 6, 9, 21, 13, 41, 34, 3, 28, 35, 37, 39, 51, 20, 27, 1, 7, 34, 39, 33, 38, 31, 30, 34, 18, 58, 5, 7, 12, 18, 21, 24, 5, 33, 31, 26, 8, 1, 7, 9, 13, 17, 23, 25, 49, 38, 11, 32, 26, 50, 5, 11, 17, 36, 1, 33, 4, 8, 11, 28, 41, 7, 21, 66, 58, 64, 51, 63, 31, 45, 50, 58, 61, 60, 63, 56, 54, 52, 46, 44, 41, 36, 35, 35, 31, 31, 59, 30, 30, 30, 27, 30, 24, 22, 30, 22, 22, 20, 15, 15, 10, 9, 49, 8, 14, 8, 8, 14, 5, 3, 1, 3, 12, 14, 17, 24, 1, 29, 3, 1, 13, 23, 7, 18, 21, 3, 52, 31, 28, 31, 45, 30, 15, 20, 25, 34, 47, 3, 17, 61, 22, 9, 25, 30, 56, 15, 24, 31, 16, 32, 37, 67, 10, 20, 22, 25, 5, 17, 47, 50, 26, 33, 13, 68, 26, 65, 18, 32, 1, 53, 66, 3, 32, 1, 1, 42, 1, 8, 51, 1, 13, 16, 6, 15, 3, 6, 9, 11, 17, 27, 1, 12, 20, 23, 6, 8, 10, 11, 17, 28, 33, 19, 26, 28, 15, 1, 5, 6, 21, 29, 1, 33, 3, 23, 3, 10, 8, 12, 5, 11, 20, 23, 25, 24, 39, 12, 14, 33, 28, 12, 55, 6, 4, 34, 10, 25, 31, 28, 10, 27, 3, 28, 11, 61, 21, 21, 16, 9, 13, 18, 3, 13, 1, 14, 3, 1, 1, 8, 32, 24, 37, 4, 37, 34, 25, 56, 10, 15, 67, 33, 18, 35, 10, 35, 48, 16, 38, 12, 35, 28, 14, 1, 1, 18, 6, 8, 12, 15, 17, 21, 5, 29, 40, 11, 25, 1, 1, 9, 25, 51, 1, 14, 31, 19, 35, 11, 24, 1, 3, 5, 4, 8, 15, 59, 76, 30, 32, 54, 53, 22, 5, 14, 18, 20, 58, 64, 35, 57, 30, 17, 54, 42, 23, 35, 32, 4, 8, 61, 14, 20, 23, 21, 27, 37, 68, 8, 17, 27, 29, 14, 24, 58, 42, 68, 4, 8, 22, 27, 7, 10, 16, 28, 56, 34, 20, 28, 21, 27, 43, 3, 30, 14, 23, 64, 29, 33, 30, 38, 13, 3, 55, 2, 1, 1, 47, 52, 26, 6, 18, 16, 18, 8, 1, 23, 48, 36, 51, 49, 41, 61, 16, 29, 11, 76, 6, 24, 19, 25, 63, 25, 35, 3, 25, 7, 25, 27, 4, 12, 26, 29, 26, 49, 5, 23, 8, 30, 70, 16, 24, 25, 66, 1, 11, 26, 31, 29, 34, 1, 29, 4, 2, 35, 63, 69, 1, 16, 10, 13, 30, 43, 7, 15, 21, 26, 62, 34, 14, 11, 14, 27, 14, 5, 13, 37, 39, 41, 32, 32, 5, 8, 28, 23, 20, 2, 3, 15, 24, 15, 15, 17, 58, 56, 44, 32, 20, 67, 26, 19, 25, 9, 14, 19, 31, 53, 3, 39, 24, 29, 8, 18, 26, 30, 32, 31, 34, 55, 13, 23, 15, 68, 15, 14, 12, 35, 16, 49, 72, 9, 13, 25, 68, 68, 7, 22, 13, 7, 17, 26, 34, 21, 10, 12, 16, 12, 41, 43, 3, 17, 1, 1, 15, 15, 13, 22, 1, 35, 3, 7, 13, 1, 8, 35, 5, 9, 68, 58, 20, 25, 32, 4, 28, 61, 15, 50, 58, 11, 8, 18, 1, 6, 25, 24, 46, 11, 24, 27, 11, 27, 2, 4, 52, 31, 48, 23, 22, 20, 7, 61, 35, 41, 33, 5, 40, 11, 28, 37, 40, 50, 35, 50, 3, 19, 22, 21, 5, 56, 15, 27, 16, 5, 50, 46, 29, 31, 1, 61, 14, 28, 23, 19, 9, 10, 13, 17, 22, 27, 66, 3, 12, 42, 11, 1, 22, 19, 10, 26, 33, 4, 11, 15, 17, 25, 26, 57, 26, 37, 5, 20, 24, 4, 10, 38, 31, 18, 28, 1, 8, 10, 13, 23, 59, 56, 56, 13, 50, 60, 3, 46, 18, 44, 62, 44, 44, 41, 48, 39, 36, 11, 29, 26, 25, 20, 16, 48, 11, 18, 11, 41, 9, 39, 7, 39, 42, 45, 53, 5, 12, 3, 49, 53, 44, 20, 4, 15, 36, 66, 15, 41, 13, 24, 1, 13, 24, 3, 9, 16, 26, 2, 5, 10, 54, 27, 66, 30, 11, 22, 2, 1, 11, 17, 57, 59, 61, 7, 8, 18, 1, 62, 13, 1, 47, 33, 28, 34, 67, 44, 72, 36, 64, 22, 50, 60, 10, 10, 7, 3, 46, 59, 11, 38, 16, 15, 13, 17, 29, 55, 22, 6, 64, 6, 12, 1, 37, 1, 36, 6, 9, 23, 5, 3, 17, 11, 22, 39, 1, 61, 26, 6, 53, 13, 24, 12, 60, 58, 53, 41, 12, 30, 6, 16, 68, 68, 5, 12, 19, 2, 3, 14, 8, 18, 2, 9, 3, 8, 17, 62, 54, 21, 28, 9, 18, 31, 1, 1, 31, 41, 20, 41, 59, 4, 15, 9, 1, 20, 30, 45, 11, 31, 33, 10, 12, 14, 16, 29, 41, 70, 20, 30, 67, 63, 61, 60, 60, 55, 51, 46, 58, 43, 42, 40, 30, 5, 29, 28, 25, 42, 17, 17, 22, 16, 11, 8, 4, 27, 4, 2, 68, 6, 15, 33, 50, 30, 17, 32, 20, 31, 33, 11, 22, 54, 56, 9, 18, 46, 56, 53, 57, 41, 1, 15, 17, 32, 4, 25, 18, 19, 21, 30, 19, 49, 9, 17, 3, 15, 3, 27, 29, 24, 33, 7, 40, 45, 18, 28, 30, 5, 11, 13, 33, 30, 3, 13, 51, 54, 2, 34, 8, 3, 11, 21, 23, 16, 12, 1, 18, 20, 7, 40, 47, 38, 45, 6, 41, 28, 17, 19, 26, 26, 3, 4, 7, 12, 5, 3, 6, 22, 3, 1, 1, 5, 9, 12, 23, 28, 3, 23, 75, 48, 47, 10, 26, 5, 13, 57, 4, 7, 10, 1, 33, 56, 33, 3, 35, 16, 34, 15, 30, 65, 5, 19, 27, 52, 31, 37, 49, 7, 13, 41, 21, 33, 10, 68, 64, 16, 20, 22, 24, 8, 20, 24, 27, 18, 28, 53, 9, 4, 26, 63, 52, 68, 59, 28, 53, 31, 9, 50, 15, 26, 32, 4, 3, 11, 17, 1, 10, 8, 66, 1, 27, 1, 56, 53, 39, 46, 34, 11, 43, 41, 18, 11, 4, 17, 27, 13, 52, 71, 4, 20, 29, 33, 22, 65, 64, 48, 2, 35, 28, 19, 30, 33, 24, 27, 15, 23, 17, 20, 23, 26, 57, 26, 32, 56, 33, 68, 12, 67, 58, 6, 9, 26, 53, 51, 41, 48, 46, 46, 2, 42, 5, 12, 31, 30, 29, 26, 25, 23, 51, 17, 50, 5, 12, 43, 9, 8, 6, 32, 5, 39, 4, 30, 2, 4, 3, 32, 3, 11, 3, 10, 20, 27, 30, 35, 38}};
        for (int[] list : lists) {

            IntArrayList inputList = IntArrayList.wrap(list);
            //System.out.println("Testing "+inputList);
            IntArrayList output = handler.deltaModTransform(inputList);
            IntList backTransformed = new IntArrayList();
            handler.decodeDeltaModTransform(output, backTransformed);
            compareLists(inputList, backTransformed);
            assertEquals("the decoded list must match the input list for " + inputList, inputList, backTransformed);
        }
    }

    private void compareLists(IntArrayList inputList, IntList backTransformed) {
        if (inputList.size()!=backTransformed.size()) {
            System.out.printf("list size differ: inputList.size=%d  backTransformed.size=%d %n",inputList.size(), backTransformed.size());
        }
        for (int i=0;i<inputList.size();i++) {
           if (inputList.getInt(i)!=backTransformed.getInt(i)) {
            System.out.printf("element differ at index=%d inputList.value=%d backTransformed.value=%d %n",   i,
                    inputList.getInt(i), backTransformed.getInt(i));
        }
        }
    }

    @Test
    public void roundTripMoreWithQualScores() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        codec.setHandler(new AlignmentCollectionHandler());
        Alignments.AlignmentCollection.Builder collection = loadCollection("test-data/seq-var-test/kevin-synth/sorted-seq-var-reads-gsnap.entries");
        addQualScores(collection);
        addToQuals(collection);
        assertRoundTripMatchExpected(codec, collection);
    }

    private void clearReadQualityScores(Alignments.AlignmentCollection.Builder collection) {
        for (int i = 0; i < collection.getAlignmentEntriesCount(); i++) {
            Alignments.AlignmentEntry.Builder element = collection.getAlignmentEntriesBuilder(i);
            element.clearReadQualityScores();
        }
    }

    private void addToQuals(Alignments.AlignmentCollection.Builder collection) {
        for (int i = 0; i < collection.getAlignmentEntriesCount(); i++) {
            Alignments.AlignmentEntry.Builder element = collection.getAlignmentEntriesBuilder(i);
            if (element.hasReadQualityScores()) {
                ByteString qualScores = element.getReadQualityScores();
                for (int j = 0; j < element.getSequenceVariationsCount(); j++) {
                    Alignments.SequenceVariation.Builder seqVarBuilder = element.getSequenceVariationsBuilder(j);
                    final String toBases = seqVarBuilder.getTo();
                    int indelOffset = 0;
                    final byte[] toQuals = new byte[toBases.length()];
                    for (int k = 0; k < toQuals.length; k++) {
                        final boolean ignoreBase = toBases.charAt(k) == '-';
                        toQuals[k] = ignoreBase ? 0 : qualScores.byteAt(seqVarBuilder.getReadIndex() - 1 + k - indelOffset);

                        if (ignoreBase) {
                            indelOffset++;
                        }
                    }
                    seqVarBuilder.setToQuality(ByteString.copyFrom(toQuals));
                }
            }

        }
    }

    private void addQualScores(Alignments.AlignmentCollection.Builder collection) {
        byte[] quals = new byte[]{1, 2, 2, 2, 3, 4, 8, 3, 2, 2, 2, 2, 9, 1, 2, 2, 2, 3, 4, 8, 3, 2, 2, 2, 2, 9, 1, 2, 2, 2, 3, 4, 8, 3, 2, 2, 2, 2, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1};
        for (int i = 0; i < collection.getAlignmentEntriesCount(); i++) {
            Alignments.AlignmentEntry.Builder element = collection.getAlignmentEntriesBuilder(i);
            element.setReadQualityScores(ByteString.copyFrom(quals));

        }

    }

    @Test
    // will not run on server.
    public void roundTripLarge() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        final AlignmentCollectionHandler handler = new AlignmentCollectionHandler();
        handler.setDebugLevel(2);
        codec.setHandler(handler);
        Alignments.AlignmentCollection.Builder collection = loadCollection("/data/rrbs/EMNWFIL.entries", 0, 100);

        assertRoundTripMatchExpected(codec, collection, false);
    }

    @Test
    // will not run on server.
    public void roundTripLarge2() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        final AlignmentCollectionHandler alignmentCollectionHandler = new AlignmentCollectionHandler();
        alignmentCollectionHandler.setEnableDomainOptimizations(true);
        codec.setHandler(alignmentCollectionHandler);
        Alignments.AlignmentCollection.Builder collection = loadCollection("/data/BAM-redo/ZHUUJKS-all-quals.entries", 0, 100000);

        assertRoundTripMatchExpected(codec, collection, false);
    }
         @Test
    // will not run on server.
    public void roundTripLargeHZ() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();
        final AlignmentCollectionHandler alignmentCollectionHandler = new AlignmentCollectionHandler();
        alignmentCollectionHandler.setEnableDomainOptimizations(true);

        codec.setHandler(alignmentCollectionHandler);
        Alignments.AlignmentCollection.Builder collection = loadCollection("/data/bug13/HZFWPTI.entries", 0, 20000);

        assertRoundTripMatchExpected(codec, collection, false);
    }
    @Test
    //roundTripExamplePairedEndDomainOptimizations
        public void roundTripExample() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();

        final AlignmentCollectionHandler alignmentCollectionHandler = new AlignmentCollectionHandler();
        alignmentCollectionHandler.setEnableDomainOptimizations(true);
        codec.setHandler(alignmentCollectionHandler);
        final Alignments.AlignmentCollection.Builder collection = loadCollection("test-data/bam/Example.entries", 0, 100000);

        assertRoundTripMatchExpected(codec, collection);
    }

    @Test
    // this test will not run on the server.
    public void roundTripUANMNXR() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();

        final AlignmentCollectionHandler alignmentCollectionHandler = new AlignmentCollectionHandler();
        alignmentCollectionHandler.setEnableDomainOptimizations(true);
        alignmentCollectionHandler.setDebugLevel(2);
        codec.setHandler(alignmentCollectionHandler);
        final Alignments.AlignmentCollection.Builder collection = loadCollection("/data/bug13/UANMNXR.entries", 0, 100);

        assertRoundTripMatchExpected(codec, collection);
    }

    @Test
    public void roundTripExamplePairedEnd() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();

        final AlignmentCollectionHandler alignmentCollectionHandler = new AlignmentCollectionHandler();
        alignmentCollectionHandler.setEnableDomainOptimizations(false);
        codec.setHandler(alignmentCollectionHandler);
        final Alignments.AlignmentCollection.Builder collection = loadCollection("test-data/bam/Example.entries", 0, 1000);

        assertRoundTripMatchExpected(codec, collection);
    }

    @Test
    public void roundTripExamplePairedSpliced() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();

        final AlignmentCollectionHandler alignmentCollectionHandler = new AlignmentCollectionHandler();
        alignmentCollectionHandler.setEnableDomainOptimizations(false);
        codec.setHandler(alignmentCollectionHandler);
        final Alignments.AlignmentCollection.Builder collection = loadCollectionNoPerm("test-data/alignment-hybrid-codec/EJOYQAZ-small.entries", 0, 1000);

        assertRoundTripMatchExpected(codec, collection);
    }

    @Test
    public void roundTripExamplePairedSplicedDomainOptimizations() throws IOException {
        final HybridChunkCodec1 codec = new HybridChunkCodec1();

        final AlignmentCollectionHandler alignmentCollectionHandler = new AlignmentCollectionHandler();
        alignmentCollectionHandler.setEnableDomainOptimizations(true);
        codec.setHandler(alignmentCollectionHandler);
        final Alignments.AlignmentCollection.Builder collection = loadCollectionNoPerm("test-data/alignment-hybrid-codec/EJOYQAZ-small.entries", 0, 1000);

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

    // @Test
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

    private Alignments.AlignmentCollection.Builder loadCollectionNoPerm(String filename, int firstElementToLoad, int maxElementsToLoad) throws IOException {
        final Alignments.AlignmentCollection.Builder collectionBuilder = Alignments.AlignmentCollection.newBuilder();
        AlignmentReaderImpl reader = new AlignmentReaderImpl(filename);
        try {
            int counter = 0;
            for (Alignments.AlignmentEntry entry : reader) {
                if (counter >= firstElementToLoad) {
                    collectionBuilder.addAlignmentEntries(entry);
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

    private void testRoundTripWithBuiltEntries(HybridChunkCodec1 codec,
                                               ObjectArrayList<Alignments.AlignmentEntry.Builder> expectedCollection,
                                               Alignments.AlignmentCollection.Builder collection,
                                               boolean addReadQual,
                                               SoftClip[] exampleClips) throws IOException {
        Alignments.AlignmentCollection.Builder expected = Alignments.AlignmentCollection.newBuilder();

        for (final Alignments.AlignmentEntry.Builder entryBuilder : expectedCollection) {

            expected.addAlignmentEntries(entryBuilder);
        }

        if (addReadQual) {

            addToQuals(expected);
            addToQuals(collection);
        } else {
            clearReadQualityScores(collection);
        }
        if (exampleClips != null) {

            addSoftClips(exampleClips, expected);
            addSoftClips(exampleClips, collection);
        }
        final ByteArrayOutputStream encoded = codec.encode(collection.build());
        Alignments.AlignmentCollection decodedCollection = (Alignments.AlignmentCollection) codec.decode(encoded.toByteArray());
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

    private Alignments.AlignmentCollection.Builder buildCollection(AlignmentExample[] builtEntries, boolean addReadQual) {
        final Alignments.AlignmentCollection.Builder collectionBuilder = Alignments.AlignmentCollection.newBuilder();
        for (Alignments.AlignmentEntry.Builder entry : buildEntriesCollection(builtEntries)) {
            addReadQuals(addReadQual, entry);
            collectionBuilder.addAlignmentEntries(entry);
        }
        return collectionBuilder;
    }

    private void addReadQuals(boolean addReadQual, Alignments.AlignmentEntry.Builder entry) {
        if (addReadQual) {

            final byte[] quals = {0x1, 0x2, 0x3, 0x42, 0x1, 0x23, 0x2, 0x3, 0x1F, 0x3, 0x2, 0x1, 0x2, 0x3, 0x42, 0x1,
                    0x23, 0x2, 0x3, 0x1F, 0x3, 0x2, 0x1, 0x2, 0x3, 0x42, 0x1, 0x23, 0x2, 0x3, 0x1F, 0x3, 0x2,};
            final ByteString scores = ByteString.copyFrom(quals);
            entry.setReadQualityScores(scores);
        }
    }

    private Alignments.AlignmentCollection.Builder buildCollectionWithQualScores(AlignmentExample[] builtEntries) {
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

    //@Test
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
        private int insertSize;


        AlignmentExample(int targetIndex, int position, int mappingQuality, int query_index, String var1_to, String var1_from, int var1_position, int var1_readIndex,
                         String var2_to, String var2_from, int var2_position, int var2_readIndex, int insertSize) {

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
            this.insertSize = insertSize;
        }
    }

    AlignmentExample[] examples = new AlignmentExample[]{
            new AlignmentExample(0, 1, 33, 1, "TA", "GG", 27, 31, "--", "TA", 3, 10, 0),
            new AlignmentExample(0, 1, 9, 1, "..", "AA", 42, 24, "T-", "G.", 5, 11, 2),
            new AlignmentExample(1, 1, 9, 1, "AA--AA", "CCCC", 13, 24, "T-", "G.", 5, 3, 34)
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
            alignmentBuilder.setQueryLength(1003);
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
            if (entry.insertSize != 0) {
                alignmentBuilder.setInsertSize(entry.insertSize);
            }
            list.add(alignmentBuilder);
        }
        return list;
    }


}
