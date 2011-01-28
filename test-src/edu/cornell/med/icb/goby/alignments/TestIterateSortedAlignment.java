/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.Test;
import org.junit.BeforeClass;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.io.File;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.ints.*;

/**
 * @author Fabien Campagne
 *         Date: Sep 3, 2010
 *         Time: 5:31:35 PM
 */
public class TestIterateSortedAlignment {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestSkipTo.class);
    private static final String BASE_TEST_DIR = "test-results/alignments-iterate-sorted";

    @BeforeClass
    public static void initializeTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Creating base test directory: " + BASE_TEST_DIR);
        }
        FileUtils.forceMkdir(new File(BASE_TEST_DIR));
    }

    @AfterClass
    public static void cleanupTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Deleting base test directory: " + BASE_TEST_DIR);
        }
        FileUtils.forceDeleteOnExit(new File(BASE_TEST_DIR));
    }

    @Test
    public void testIterateSorted() throws IOException {


        final String basename = "align-skip-to-1-concat";
        final String basenamePath = FilenameUtils.concat(BASE_TEST_DIR, basename);
        final AlignmentWriter writer =
                new AlignmentWriter(basenamePath);
        writer.setNumAlignmentEntriesPerChunk(1);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);
        Alignments.AlignmentEntry.Builder newEntry;

        newEntry = prepareAlignmentEntry(0, 1, 1, 30, false, new int[0]);
        writer.appendEntry(newEntry.build());

        newEntry = prepareAlignmentEntry(0, 1, 130, 30, false, new int[]{11, 12});
        writer.appendEntry(newEntry.build());

        newEntry = prepareAlignmentEntry(0, 1, 135, 35, false, new int[]{6});
        writer.appendEntry(newEntry.build());

        newEntry = prepareAlignmentEntry(0, 2, 1230, 30, false, new int[0]);
        writer.appendEntry(newEntry.build());


        newEntry = prepareAlignmentEntry(0, 2, 3000, 30, false, new int[0]);
        writer.appendEntry(newEntry.build());


        writer.close();
        writer.printStats(System.out);

        final Int2IntMap positionMap = new Int2IntOpenHashMap();

        IterateSortedAlignmentsListImpl iterator = new IterateSortedAlignmentsListImpl() {


            @Override
            public void processPositions(int referenceIndex, int intermediatePosition, ObjectArrayList<PositionBaseInfo> positionBaseInfos) {
                positionMap.put(intermediatePosition, positionBaseInfos.size());
                System.out.printf("position: %d listSize: %d%n", referenceIndex, positionBaseInfos.size());
            }
        };
        iterator.iterate(basenamePath);

        for (int i = 1; i < 35; i++) {
            assertEquals("position " + i, 1, positionMap.get(i));
        }
        for (int i = 130; i < 135; i++) {
            assertEquals("position " + i, 1, positionMap.get(i));
        }
        for (int i = 135; i < 165; i++) {
            assertEquals("position " + i, 2, positionMap.get(i));
        }
        for (int i = 165; i < 170; i++) {
            assertEquals("position " + i, 1, positionMap.get(i));
        }
        for (int i = 1230; i < 1265; i++) {
            assertEquals("position " + i, 1, positionMap.get(i));
        }
        for (int i = 3000; i < 3035; i++) {
            assertEquals("position " + i, 1, positionMap.get(i));
        }
    }

    @Test
    public void testIterateSortedTwoMutations() throws IOException {


        final String basename = "align-skip-to-1-concat-two-mutations";
        final String basenamePath = FilenameUtils.concat(BASE_TEST_DIR, basename);
        final AlignmentWriter writer =
                new AlignmentWriter(basenamePath);
        writer.setNumAlignmentEntriesPerChunk(1);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);
        Alignments.AlignmentEntry.Builder newEntry;


        newEntry = prepareAlignmentEntry(0, 1, 100, 3, true, new int[]{2, 31, 34}, 35);
        writer.appendEntry(newEntry.build());

        writer.close();
        writer.printStats(System.out);

        final IntSet variantReadIndices = new IntOpenHashSet();
        final IntSet variantPositionOnRef = new IntOpenHashSet();
        IterateSortedAlignmentsListImpl iterator = new IterateSortedAlignmentsListImpl() {

            public void processPositions(int referenceIndex, int intermediatePosition, ObjectArrayList<PositionBaseInfo> positionBaseInfos) {

            }

            @Override
            public void observeVariantBase(ConcatSortedAlignmentReader sortedReaders,
                                           Alignments.AlignmentEntry alignmentEntry, Int2ObjectMap<ObjectArrayList<PositionBaseInfo>> positionToBases,
                                           Alignments.SequenceVariation var, char toChar, char fromChar,
                                           int currentReferenceIndex, int currentRefPosition, int currentReadIndex) {
                variantReadIndices.add(currentReadIndex);
                variantPositionOnRef.add(currentRefPosition);
            }


        };
        iterator.iterate(basenamePath);

        assertTrue(variantReadIndices.contains(34));
        assertTrue(variantReadIndices.contains(2));
        assertTrue(variantReadIndices.contains(5));
        assertTrue(variantPositionOnRef.contains(101));
        assertTrue(variantPositionOnRef.contains(130));
        assertTrue(variantPositionOnRef.contains(133));

    }

    @Test
    public void testIterateSortedTwoMutationsSmall() throws IOException {


        final String basename = "align-skip-to-1-concat-two-mutations-small";
        final String basenamePath = FilenameUtils.concat(BASE_TEST_DIR, basename);
        final AlignmentWriter writer =
                new AlignmentWriter(basenamePath);
        writer.setNumAlignmentEntriesPerChunk(1);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);
        Alignments.AlignmentEntry.Builder newEntry;


        newEntry = prepareAlignmentEntry(0, 1, 100, 3, true, new int[]{2, 4}, 5);
        writer.appendEntry(newEntry.build());

        writer.close();
        writer.printStats(System.out);

        final IntSet variantReadIndices = new IntOpenHashSet();
        final IntSet variantPositionOnRef = new IntOpenHashSet();
        IterateSortedAlignmentsListImpl iterator = new IterateSortedAlignmentsListImpl() {

            public void processPositions(int referenceIndex, int intermediatePosition, ObjectArrayList<PositionBaseInfo> positionBaseInfos) {

            }

            @Override
            public void observeVariantBase(ConcatSortedAlignmentReader sortedReaders,
                                           Alignments.AlignmentEntry alignmentEntry, Int2ObjectMap<ObjectArrayList<PositionBaseInfo>> positionToBases,
                                           Alignments.SequenceVariation var, char toChar, char fromChar,
                                           int currentReferenceIndex, int currentRefPosition, int currentReadIndex) {
                variantReadIndices.add(currentReadIndex);
                variantPositionOnRef.add(currentRefPosition);
            }


        };
        iterator.iterate(basenamePath);


        assertTrue(variantReadIndices.contains(2));
        assertTrue(variantReadIndices.contains(4));
        assertTrue(variantPositionOnRef.contains(101));
        assertTrue(variantPositionOnRef.contains(103));


    }

    private Alignments.AlignmentEntry.Builder prepareAlignmentEntry(final int queryIndex, final int targetIndex,
                                                                    final int position,
                                                                    final float score, final boolean matchesReverseStrand,
                                                                    int[] variationIndices) {
        return prepareAlignmentEntry(queryIndex, targetIndex, position, score, matchesReverseStrand, variationIndices, 35);
    }

    private Alignments.AlignmentEntry.Builder prepareAlignmentEntry(final int queryIndex, final int targetIndex,
                                                                    final int position,
                                                                    final float score, final boolean matchesReverseStrand,
                                                                    int[] variationIndices, int queryLength) {
        Alignments.AlignmentEntry.Builder newEntry = Alignments.AlignmentEntry.newBuilder();
        newEntry.setQueryIndex(queryIndex);
        newEntry.setTargetIndex(targetIndex);
        newEntry.setScore(score);
        newEntry.setPosition(position);
        newEntry.setMatchingReverseStrand(matchesReverseStrand);
        newEntry.setMultiplicity(1);
        newEntry.setQueryLength(queryLength);

        for (int variaIndex : variationIndices) {
            Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder();
            varBuilder.setFrom("A");
            varBuilder.setTo("C");
            varBuilder.setReadIndex(1);
            varBuilder.setPosition(variaIndex);

            newEntry.addSequenceVariations(varBuilder.build());
        }
        return newEntry;
    }

    private Alignments.AlignmentEntry.Builder prepareAlignmentEntryWithReferenceInsertion(final int queryIndex,
                                                                                          final int targetIndex,
                                                                                          final int position,
                                                                                          final float score,
                                                                                          final boolean matchesReverseStrand,
                                                                                          int[] variationIndices) {

        Alignments.AlignmentEntry.Builder newEntry = Alignments.AlignmentEntry.newBuilder();
        newEntry.setQueryIndex(queryIndex);
        newEntry.setTargetIndex(targetIndex);
        newEntry.setScore(score);
        newEntry.setPosition(position);
        newEntry.setMatchingReverseStrand(matchesReverseStrand);
        newEntry.setMultiplicity(1);
        newEntry.setQueryLength(35);


        for (int i = 0; i < variationIndices.length / 2; i += 2) {
            Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder();
            StringBuffer insertion = new StringBuffer();
            StringBuffer gaps = new StringBuffer();
            final int gapLength = variationIndices[i + 1];
            for (int j = 0; j < gapLength; j++) {
                gaps.append('-');
                insertion.append('C');
            }
            varBuilder.setTo(gaps.toString());
            varBuilder.setFrom(insertion.toString());
            final int variationIndex = variationIndices[i + 0];
            varBuilder.setPosition(variationIndex);
            varBuilder.setReadIndex(variationIndex);


            newEntry.addSequenceVariations(varBuilder.build());
        }
        return newEntry;
    }

    private Alignments.AlignmentEntry.Builder prepareAlignmentEntryWithReadInsertion(final int queryIndex, final int targetIndex,
                                                                                     final int position,
                                                                                     final float score, final boolean matchesReverseStrand,
                                                                                     int[] variationIndices) {
        Alignments.AlignmentEntry.Builder newEntry = Alignments.AlignmentEntry.newBuilder();
        newEntry.setQueryIndex(queryIndex);
        newEntry.setTargetIndex(targetIndex);
        newEntry.setScore(score);
        newEntry.setPosition(position);
        newEntry.setMatchingReverseStrand(matchesReverseStrand);
        newEntry.setMultiplicity(1);
        newEntry.setQueryLength(35);


        for (int i = 0; i < variationIndices.length / 2; i += 2) {
            Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder();
            StringBuffer insertion = new StringBuffer();
            StringBuffer gaps = new StringBuffer();
            final int gapLength = variationIndices[i + 1];
            for (int j = 0; j < gapLength; j++) {
                gaps.append('-');
                insertion.append('C');
            }
            varBuilder.setFrom(gaps.toString());
            varBuilder.setTo(insertion.toString());
            final int variationIndex = variationIndices[i + 0];
            varBuilder.setPosition(variationIndex);
            varBuilder.setReadIndex(variationIndex);


            newEntry.addSequenceVariations(varBuilder.build());
        }
        return newEntry;
    }

    @Test
    public void testIterateSortedReferenceInsertions() throws IOException {


        final String basename = "align-skip-to-1-concat";
        final String basenamePath = FilenameUtils.concat(BASE_TEST_DIR, basename);
        final AlignmentWriter writer =
                new AlignmentWriter(basenamePath);
        writer.setNumAlignmentEntriesPerChunk(1);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);
        Alignments.AlignmentEntry.Builder newEntry;

        newEntry = prepareAlignmentEntryWithReferenceInsertion(0, 1, 100, 30, false, new int[]{4, 3});
        writer.appendEntry(newEntry.build());


        writer.close();
        writer.printStats(System.out);

        final Int2IntMap positionMap = new Int2IntOpenHashMap();
        final IntSet variantReadIndices = new IntOpenHashSet();

        IterateSortedAlignmentsListImpl iterator = new IterateSortedAlignmentsListImpl() {
            @Override
            public void observeVariantBase(ConcatSortedAlignmentReader sortedReaders, Alignments.AlignmentEntry alignmentEntry, Int2ObjectMap<ObjectArrayList<PositionBaseInfo>> positionToBases, Alignments.SequenceVariation var, char toChar, char fromChar, int currentReferenceIndex, int currentRefPosition, int currentReadIndex) {
                variantReadIndices.add(currentReadIndex);
            }

            @Override
            public void processPositions(int referenceIndex, int intermediatePosition, ObjectArrayList<PositionBaseInfo> positionBaseInfos) {
                int coverage = 0;
                for (PositionBaseInfo info : positionBaseInfos) {
                    coverage += info.to != '-' ? 1 : 0;
                }
                positionMap.put(intermediatePosition, coverage);
                System.out.printf("position: %d listSize: %d%n", referenceIndex, coverage);
            }
        };
        iterator.iterate(basenamePath);

        for (int i = 0; i < 100; i++) {
            assertEquals("position " + i, 0, positionMap.get(i));
        }
        for (int i = 100; i < 103; i++) {
            assertEquals("position " + i, 1, positionMap.get(i));
        }
        for (int i = 103; i < 106; i++) {
            assertEquals("position " + i, 0, positionMap.get(i));
        }
        for (int i = 106; i < 137; i++) {
            assertEquals("position " + i, 1, positionMap.get(i));
        }
        for (int i = 138; i < 150; i++) {
            assertEquals("position " + i, 0, positionMap.get(i));
        }

    }

    @Test
    public void testIterateSortedWithReadInsertions() throws IOException {


        final String basename = "align-skip-to-1-concat";
        final String basenamePath = FilenameUtils.concat(BASE_TEST_DIR, basename);
        final AlignmentWriter writer =
                new AlignmentWriter(basenamePath);
        writer.setNumAlignmentEntriesPerChunk(1);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);
        Alignments.AlignmentEntry.Builder newEntry;

        newEntry = prepareAlignmentEntryWithReadInsertion(0, 1, 100, 30, false, new int[]{4, 3});
        writer.appendEntry(newEntry.build());


        writer.close();
        writer.printStats(System.out);

        final Int2IntMap positionMap = new Int2IntOpenHashMap();
        final Int2IntMap readIndexMap = new Int2IntOpenHashMap();

        IterateSortedAlignmentsListImpl iterator = new IterateSortedAlignmentsListImpl() {


            @Override
            public void processPositions(int referenceIndex, int intermediatePosition, ObjectArrayList<PositionBaseInfo> positionBaseInfos) {
                int coverage = 0;
                for (PositionBaseInfo info : positionBaseInfos) {
                    coverage += info.from != '-' ? 1 : 0;
                    readIndexMap.put(info.readIndex, readIndexMap.get(info.readIndex) + 1);
                }
                positionMap.put(intermediatePosition, coverage);
                //    System.out.printf("position: %d listSize: %d%n", position, coverage);
            }
        };
        iterator.iterate(basenamePath);
        // check coverage on ref:
        for (int i = 0; i < 100; i++) {
            assertEquals("position " + i, 0, positionMap.get(i));
        }
        for (int i = 100; i < 132; i++) {

            assertEquals("position " + i, 1, positionMap.get(i));
        }
        //check read index:
        for (int i = 1; i < 35; i++) {
            assertEquals("read-index " + i, 1, readIndexMap.get(i));

        }

    }


    @Test
    public void testIterateSortedTwoTargetSequences() throws IOException {


        final String basename = "align-skip-to-1-contact-two-targets";
        final String basenamePath = FilenameUtils.concat(BASE_TEST_DIR, basename);
        final AlignmentWriter writer =
                new AlignmentWriter(basenamePath);
        writer.setNumAlignmentEntriesPerChunk(1);

        final int numTargets = 3;
        final int[] targetLengths = new int[numTargets];

        for (int referenceIndex = 0; referenceIndex < numTargets; referenceIndex++) {
            targetLengths[referenceIndex] = 1000;
        }
        writer.setTargetLengths(targetLengths);
        // we write this alignment sorted:

        writer.setSorted(true);
        Alignments.AlignmentEntry.Builder newEntry;

        newEntry = prepareAlignmentEntry(0, 1, 1, 30, false, new int[0]);
        writer.appendEntry(newEntry.build());

        newEntry = prepareAlignmentEntry(0, 1, 130, 30, false, new int[]{11, 12});
        writer.appendEntry(newEntry.build());

        newEntry = prepareAlignmentEntry(0, 1, 135, 35, false, new int[]{6});
        writer.appendEntry(newEntry.build());

        newEntry = prepareAlignmentEntry(0, 2, 1, 30, false, new int[0]);
        writer.appendEntry(newEntry.build());


        newEntry = prepareAlignmentEntry(0, 2, 130, 30, false, new int[0]);
        writer.appendEntry(newEntry.build());


        writer.close();
        writer.printStats(System.out);

        final Int2IntMap positionMap = new Int2IntOpenHashMap();

        IterateSortedAlignmentsListImpl iterator = new IterateSortedAlignmentsListImpl() {


            @Override
            public void processPositions(int referenceIndex, int intermediatePosition, ObjectArrayList<PositionBaseInfo> positionBaseInfos) {
                if (referenceIndex == 1) {
             // record only reference 1 matches. 
                    positionMap.put(intermediatePosition, positionBaseInfos.size());
                }
                System.out.printf("position: %d listSize: %d%n", referenceIndex, positionBaseInfos.size());

            }
        };
        iterator.iterate(basenamePath);

        for (int i = 1; i < 35; i++) {
            assertEquals("position " + i, 1, positionMap.get(i));
        }
        for (int i = 130; i < 135; i++) {
            assertEquals("position " + i, 1, positionMap.get(i));
        }
        for (int i = 135; i < 165; i++) {
            assertEquals("position " + i, 2, positionMap.get(i));
        }
        for (int i = 165; i < 170; i++) {
            assertEquals("position " + i, 1, positionMap.get(i));
        }

    }

}
