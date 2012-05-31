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
import edu.cornell.med.icb.goby.alignments.processors.InfoForTarget;
import edu.cornell.med.icb.goby.alignments.processors.ObservedIndel;
import edu.cornell.med.icb.goby.alignments.processors.RealignmentProcessor;
import edu.cornell.med.icb.goby.modes.AbstractAlignmentToCompactMode;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectListIterator;
import it.unimi.dsi.lang.MutableString;

import static org.junit.Assert.*;

import org.junit.Test;

import javax.swing.*;
import java.io.IOException;
import java.util.Arrays;

/**
 * Test on the fly realignment around indels. See description of test-cases at
 * https://docs.google.com/document/d/1hjptov0YuMiulG3Gd-0pLrQz_pvRXozzYwJX7Y8-UJY/edit?hl=en#
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2011
 *         Time: 12:04:07 PM
 */
public class TestRealignmentProcessor {

    private ObjectListIterator<Alignments.AlignmentEntry> buildList1() {
        ObjectList<Alignments.AlignmentEntry> list = new ObjectArrayList<Alignments.AlignmentEntry>();
        addIndel(list, 10, 50, 3, 20); // indel in middle of the read, would be correctly aligned.
        addMutation(list, 1, 50, 20, 'A', 'T');
        return list.iterator();
    }

    private String[] list2Refs() {
        return new String[]{"ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA"};
    }


    private ObjectListIterator<Alignments.AlignmentEntry> buildList2() {
        ObjectList<Alignments.AlignmentEntry> list = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "ACTGACTGACTGAATTACTA");  // this read should be realigned to the right
        addEntry(list, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAA----TTACTAG"); // this read carries the candidate indel
        /*
queryIndex=0:
var.readIndex=9
var.position=10
var.from="CTAG"
var.to="----"
entry.position=5

queryIndex=1:
var.readIndex=14 (defined as the position before the left-most '-')
var.position=15 (defined as the position of the left-most '-')
var.from="CTAG"
var.to="----"
entry.position=0
         */
        return list.iterator();
    }

    private ObjectListIterator<Alignments.AlignmentEntry> buildList2Left() {
        ObjectList<Alignments.AlignmentEntry> list = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "ACTGACTGACTGAACTA");  // this read should be realigned to the right
        addEntry(list, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "ACT------GACTGACTGAACTA"); // this read carries the candidate indel

        return list.iterator();
    }

    private ObjectListIterator<Alignments.AlignmentEntry> buildListDifferentTargets() {
        ObjectList<Alignments.AlignmentEntry> list = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list, 0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "ACTGACTGACTGAATTACTA");  // this read should be realigned to the right
        addEntry(list, 0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAA----TTACTAG"); // this read carries the candidate indel

        addEntry(list, 1, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "ACTGACTGACTGAATTACTA");  // this read should be realigned to the right
        addEntry(list, 1, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAA----TTACTAG"); // this read carries the candidate indel

        addEntry(list, 4, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "ACTGACTGACTGAATTACTA");  // this read should be realigned to the right
        addEntry(list, 4, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAA----TTACTAG"); // this read carries the candidate indel

        return list.iterator();
    }

    @Test
    public void pushEntryToIndelState() {

        ObjectListIterator<Alignments.AlignmentEntry> list = buildList0();
        RealignmentProcessor realigner = new RealignmentProcessor(list);
        InfoForTarget tinfo = new InfoForTarget(0);
        final Alignments.AlignmentEntry entry = list.next();
        System.out.println("entry: " + entry);
        realigner.pushEntryToPool(tinfo, 0, entry);
        assertEquals(2, getFirst(tinfo.potentialIndels).getStart());
        assertEquals(5, getFirst(tinfo.potentialIndels).getEnd());
        assertEquals(3, getFirst(tinfo.potentialIndels).positionSpan());

    }

    private ObservedIndel getFirst(ObjectAVLTreeSet<ObservedIndel> potentialIndels) {
        return potentialIndels.first();
    }

    private ObjectListIterator<Alignments.AlignmentEntry> buildList0() {
        ObjectList<Alignments.AlignmentEntry> list = new ObjectArrayList<Alignments.AlignmentEntry>();
        addIndel(list, 0, 10, 3, 3); // indel at position 2+1, length 3 AA---CC  start=2 end=5
        return list.iterator();
    }

    @Test
    public void testFirstRealign() throws IOException {

        ObjectListIterator<Alignments.AlignmentEntry> list = buildList1();
        RealignmentProcessor realigner = new RealignmentProcessor(list) {
            int entryIndex = 0;

            @Override
            public void pushEntryToPool(InfoForTarget tinfo, int position, Alignments.AlignmentEntry entry) {
                super.pushEntryToPool(tinfo, position, entry);
                entryIndex++;
                if (entryIndex == 1) {
                    assertFalse(tinfo.positionsWithSpanningIndel.contains(29 - 1));
                    assertTrue(tinfo.positionsWithSpanningIndel.contains(30 - 1));
                    assertTrue(tinfo.positionsWithSpanningIndel.contains(31 - 1));
                    assertTrue(tinfo.positionsWithSpanningIndel.contains(32 - 1));
                    assertFalse(tinfo.positionsWithSpanningIndel.contains(33 - 1));

                    assertEquals(30 - 1, getFirst(tinfo.potentialIndels).getStart());
                    assertEquals(33 - 1, getFirst(tinfo.potentialIndels).getEnd());
                }
            }
        };
        Alignments.AlignmentEntry entry;
        while ((entry = realigner.nextRealignedEntry(0, 0)) != null) {
            System.out.println("entry:"
                    + entry);
        }
    }

    @Test
    public void testProcessDifferentTargets() throws IOException {
        ObjectListIterator<Alignments.AlignmentEntry> list = buildListDifferentTargets();
        RealignmentProcessor realigner = new RealignmentProcessor(list);
        Alignments.AlignmentEntry entry;
        int[] count = new int[5];
        while ((entry = realigner.nextRealignedEntry(0, 0)) != null) {
            System.out.printf("processing entry on target %d at position %d %n",
                    entry.getTargetIndex(), entry.getPosition());
            count[entry.getTargetIndex()]++;

        }
        assertEquals(2, count[0]);
        assertEquals(2, count[1]);
        assertEquals(2, count[4]);

    }

    /**
     * Test case 1
     */
    @Test
    public void testRealignToTheRight() throws IOException {


        ObjectListIterator<Alignments.AlignmentEntry> list = buildList2();
        RealignmentProcessor realigner = new RealignmentProcessor(list);
        realigner.setGenome(new RandomAccessSequenceTestSupport(list2Refs()));
        Alignments.AlignmentEntry entry;
        while ((entry = realigner.nextRealignedEntry(0, 0)) != null) {

            if (entry.getQueryIndex() == 0) {
                Alignments.SequenceVariation var = entry.getSequenceVariations(0);
                assertEquals("CTAG", var.getFrom());
                assertEquals("----", var.getTo());
                assertEquals(15, var.getPosition());
                assertEquals(15, var.getReadIndex());

            }

        }
    }

    @Test
    /**
     * Test case 2
     */
    public void testRealignToTheRightReverseStrand() throws IOException {

        ObjectList<Alignments.AlignmentEntry> list1 = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list1, 0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "ACTGACTGACTGAATTACTA", true);  // this read should be realigned to the right
        addEntry(list1, 0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAA----TTACTAG"); // this read carries the candidate indel

        ObjectListIterator<Alignments.AlignmentEntry> list = list1.iterator();

        RealignmentProcessor realigner = new RealignmentProcessor(list);
        realigner.setGenome(new RandomAccessSequenceTestSupport(list2Refs()));
        Alignments.AlignmentEntry entry;
        while ((entry = realigner.nextRealignedEntry(0, 0)) != null) {

            if (entry.getQueryIndex() == 0) {
                assertTrue("entry must match reverse strand", entry.getMatchingReverseStrand());
                Alignments.SequenceVariation var = entry.getSequenceVariations(0);
                assertEquals("CTAG", var.getFrom());
                assertEquals("----", var.getTo());
                assertEquals(15, var.getPosition());
                assertEquals(7, var.getReadIndex());

            }
        }
    }

    /**
     * Test case 3
     */
    @Test
    public void testRealignToTheLeft() throws IOException {


        ObjectList<Alignments.AlignmentEntry> list1 = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list1, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "        CCTGAACTAGTTACTAGCTA");  // this read should be realigned to the right
        addEntry(list1, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "AC-------CTGAACTAGTTACTAGCTA"); // this read carries the candidate indel
        // ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA
        // ACT------GACTGACTGAACTA
        ObjectListIterator<Alignments.AlignmentEntry> iterator = list1.iterator();
        RealignmentProcessor realigner = new RealignmentProcessor(iterator);
        realigner.setGenome(new RandomAccessSequenceTestSupport(list2Refs()));
        Alignments.AlignmentEntry entry;
        while ((entry = realigner.nextRealignedEntry(0, 0)) != null) {

            if (entry.getQueryIndex() == 0) {
                Alignments.SequenceVariation var = entry.getSequenceVariations(0);
                assertEquals("TGACTGA", var.getFrom());
                assertEquals("-------", var.getTo());
                assertEquals(2, var.getPosition());

            }
        }
    }

    @Test
    /**
     * Test case 3 & 4
     */
    public void testRealignToTheLeft2() throws IOException {


        ObjectList<Alignments.AlignmentEntry> list1 = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list1, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "AC-------CTGAACTAGTTACTAGCTA"); // this read carries the candidate indel
        addEntry(list1, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "       ACCTGAACTAGTTACTAGCTA");  // this read should be realigned to the right
        addEntry(list1, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "        CCTGAACTAGTTACTAGCTA");  // this read should be realigned to the right
        // ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA
        // ACT------GACTGACTGAACTA
        ObjectListIterator<Alignments.AlignmentEntry> iterator = list1.iterator();
        RealignmentProcessor realigner = new RealignmentProcessor(iterator);
        realigner.setGenome(new RandomAccessSequenceTestSupport(list2Refs()));
        Alignments.AlignmentEntry entry;
        while ((entry = realigner.nextRealignedEntry(0, 0)) != null) {

            if (entry.getQueryIndex() == 2) {
                Alignments.SequenceVariation var = entry.getSequenceVariations(0);
                assertEquals("TGACTGA", var.getFrom());
                assertEquals("-------", var.getTo());
                assertEquals(2, var.getPosition());
                assertEquals(2, var.getReadIndex());

            } else if (entry.getQueryIndex() == 1) {
                Alignments.SequenceVariation var = entry.getSequenceVariations(0);
                assertEquals("TGACTGA", var.getFrom());
                assertEquals("-------", var.getTo());
                assertEquals(3, var.getPosition());
                assertEquals(3, var.getReadIndex());

            }

        }
    }

    @Test
    /**
     * Test case 5
     */
    public void testRealignToTheLeftReverseStrand3() throws IOException {


        ObjectList<Alignments.AlignmentEntry> list1 = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list1, 0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "AC-------CTGAACTAGTTACTAGCTA"); // this read carries the candidate indel
        addEntry(list1, 0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "        CCTGAACTAGTTACTAGCTA", true);  // this read should be realigned to the right
        // ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA
        // ACT------GACTGACTGAACTA
        ObjectListIterator<Alignments.AlignmentEntry> iterator = list1.iterator();
        RealignmentProcessor realigner = new RealignmentProcessor(iterator);
        realigner.setGenome(new RandomAccessSequenceTestSupport(list2Refs()));
        Alignments.AlignmentEntry entry;
        while ((entry = realigner.nextRealignedEntry(0, 0)) != null) {

            if (entry.getQueryIndex() == 1) {
                assertTrue(entry.getMatchingReverseStrand());
                Alignments.SequenceVariation var = entry.getSequenceVariations(0);
                assertEquals("TGACTGA", var.getFrom());
                assertEquals("-------", var.getTo());
                assertEquals(2, var.getPosition());
                assertEquals(26, var.getReadIndex());

            }

        }
    }

    /**
     * Test case 6
     *
     * @throws java.io.IOException
     */
    @Test
    public void testCase6() throws IOException {


        ObjectList<Alignments.AlignmentEntry> list1 = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list1, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "ACTGACTGACTGAATTACTAA");  // this read should be realigned to the right
        addEntry(list1, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAA----TTACTAG"); // this read carries the candidate indel

        ObjectListIterator<Alignments.AlignmentEntry> list = list1.iterator();
        RealignmentProcessor realigner = new RealignmentProcessor(list);
        realigner.setGenome(new RandomAccessSequenceTestSupport(list2Refs()));
        Alignments.AlignmentEntry entry;
        while ((entry = realigner.nextRealignedEntry(0, 0)) != null) {

            if (entry.getQueryIndex() == 0) {
                assertEquals(2, entry.getSequenceVariationsCount());
                Alignments.SequenceVariation var1 = entry.getSequenceVariations(1);
                assertEquals("CTAG", var1.getFrom());
                assertEquals("----", var1.getTo());
                assertEquals(15, var1.getPosition());
                assertEquals(15, var1.getReadIndex());
                Alignments.SequenceVariation var2 = entry.getSequenceVariations(0);
                assertEquals("G", var2.getFrom());
                assertEquals("A", var2.getTo());
                assertEquals(25, var2.getPosition());
                assertEquals(25, var2.getReadIndex());

            }

        }
    }


    /**
     * Test case 7  TODO: enable this test for read insertion.
     */
    public void testCase7() throws IOException {


        ObjectList<Alignments.AlignmentEntry> list1 = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list1, 0, "ACTGACTGACTGAATTACTAGCTAAAGTTA", "     CTGACTGAACTAGTTACTAG");  // this read should be realigned to the right
        addEntry(list1, 0, "ACTGACTGACTGAA----TTACTAGCTAAAGTTA", "     CTGACTGAACTAGTTACTAG"); // this read carries the candidate read insertion
        // ACTGACTGACTGAA----TTACTAGCTAAAGTTA ref
        //      CTGACTGAACTAGTTACTAG          read with insertion.    
        ObjectListIterator<Alignments.AlignmentEntry> iterator = list1.iterator();
        RealignmentProcessor realigner = new RealignmentProcessor(iterator);
        realigner.setGenome(new RandomAccessSequenceTestSupport(list2Refs()));
        Alignments.AlignmentEntry entry;
        while ((entry = realigner.nextRealignedEntry(0, 0)) != null) {

            if (entry.getQueryIndex() == 0) {
                assertFalse(entry.getMatchingReverseStrand());
                Alignments.SequenceVariation var = entry.getSequenceVariations(0);
                assertEquals("----", var.getFrom());
                assertEquals("CTAG", var.getTo());
                assertEquals(10, var.getPosition());
                assertEquals(10, var.getReadIndex());

            }
            System.out.println("entry:"
                    + entry);
        }
    }

    private ObjectListIterator<Alignments.AlignmentEntry> buildList3() {
        ObjectList<Alignments.AlignmentEntry> list = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list, 0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "ACTGACTGACTGAATTACTA", true);  // this read should be realigned to the right
        addEntry(list, 0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAA----TTACTAG"); // this read carries the candidate indel

        return list.iterator();
    }

    @Test
    public void testRealignmentScore1() {
        ObjectListIterator<Alignments.AlignmentEntry> iterator = new ObjectArrayList().iterator();
        RealignmentProcessor realigner = new RealignmentProcessor(iterator);
        RandomAccessSequenceInterface genome = new RandomAccessSequenceTestSupport(list2Refs());
        Alignments.AlignmentEntry entry = makeEntry(0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAATTACTA", true).build();
        ObservedIndel indel = new ObservedIndel(14, 18, "CTAG", "----");
        assertEquals(3, realigner.score(entry, indel, true, 0, genome));

        entry = makeEntry(0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAATTACTAG", true).build();
        indel = new ObservedIndel(14, 18, "CTAG", "----");
        assertEquals(4, realigner.score(entry, indel, true, 0, genome));
    }

    @Test
    public void testRealignmentScore2() {
        ObjectListIterator<Alignments.AlignmentEntry> iterator = new ObjectArrayList().iterator();
        RealignmentProcessor realigner = new RealignmentProcessor(iterator);
        RandomAccessSequenceInterface genome = new RandomAccessSequenceTestSupport(list2Refs());
        Alignments.AlignmentEntry entry = makeEntry(0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAATTAGTA", true).build();
        ObservedIndel indel = new ObservedIndel(15, 19, "CTAG", "----");
        // expected score is 1 because despite the +3 you get for inserting the indel you have one base mismatch between the ref and the read.
        assertEquals(1, realigner.score(entry, indel, true, 0, genome));
    }

    private int queryIndex = 0;


    private void addEntry(ObjectList<Alignments.AlignmentEntry> list, int targetIndex, String reference, String read) {
        addEntry(list, targetIndex, reference, read, false);
    }

    private void addEntry(ObjectList<Alignments.AlignmentEntry> list, int targetIndex, String reference, String read, boolean reverseStrand) {
        Alignments.AlignmentEntry.Builder entry = makeEntry(targetIndex, reference, read, reverseStrand);
        list.add(entry.build());
    }

    private Alignments.AlignmentEntry.Builder makeEntry(int targetIndex, String reference, String read, boolean reverseStrand) {
        Alignments.AlignmentEntry.Builder entry = Alignments.AlignmentEntry.newBuilder();

        entry.setTargetIndex(targetIndex);
        entry.setMatchingReverseStrand(reverseStrand);
        int position = -1;
        while (read.charAt(++position) == ' ') ;
        entry.setPosition(position);
        final String querySequenceNoGaps = read.trim().
                replaceAll("-", "").
                replaceAll(" ", "");
        int queryLength = querySequenceNoGaps.length();
        final String querySequenceOnlyGaps = read.trim().
                replaceAll("A", "").
                replaceAll("C", "").
                replaceAll("T", "").replaceAll("G", "");
        int deletionLength = querySequenceOnlyGaps.length();
        int alignmentLength = queryLength + deletionLength;
        entry.setQueryLength(queryLength);
        entry.setQueryAlignedLength(alignmentLength);
        entry.setTargetAlignedLength(read.trim().length());
        entry.setQueryIndex(queryIndex++);

        /*
       final Alignments.AlignmentEntry.Builder currentEntry, final int alignmentLength,
                                                final MutableString referenceSequence,
                                                final MutableString readSequence,
                                                final int readStartPosition,
                                                final int queryLength, final boolean reverseStrand,
                                                byte[] baseQualities
        */
        final MutableString referenceMutable = new MutableString(reference.substring(position, position + alignmentLength));
        final MutableString readMutable = new MutableString(read.trim());
        AbstractAlignmentToCompactMode.extractSequenceVariations(entry, alignmentLength,
                referenceMutable,
                readMutable,
                0, queryLength, false, baseQual(read));
        return entry;
    }

    private void addEntry(ObjectList<Alignments.AlignmentEntry> list, String reference, String read) {
        addEntry(list, 0, reference, read);
    }

    private byte[] baseQual(String read) {
        byte[] result = new byte[read.length()];
        Arrays.fill(result, (byte) 40);
        return result;
    }

    @Test
    public void testIndel() {
        ObjectList<Alignments.AlignmentEntry> list = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAA----TTACTAG");
        Alignments.AlignmentEntry entry = list.get(0);
        assertEquals("----", entry.getSequenceVariations(0).getTo());
        assertEquals("CTAG", entry.getSequenceVariations(0).getFrom());
    }

    @Test
    public void testMutations() {
        ObjectList<Alignments.AlignmentEntry> list = new ObjectArrayList<Alignments.AlignmentEntry>();
        addEntry(list, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "ACTGACTGACTGAATTACTA");
        Alignments.AlignmentEntry entry = list.get(0);
        assertEquals(15, entry.getSequenceVariations(0).getPosition());
        assertEquals("C", entry.getSequenceVariations(0).getFrom());
        assertEquals("T", entry.getSequenceVariations(0).getTo());

        assertEquals(18, entry.getSequenceVariations(1).getPosition());
        assertEquals("G", entry.getSequenceVariations(1).getFrom());
        assertEquals("C", entry.getSequenceVariations(1).getTo());

        assertEquals(20, entry.getSequenceVariations(2).getPosition());
        assertEquals("T", entry.getSequenceVariations(2).getFrom());
        assertEquals("A", entry.getSequenceVariations(2).getTo());
    }

    @Test
    public void testQualityPassThrough() {
        ObjectList<Alignments.AlignmentEntry> list = new ObjectArrayList<Alignments.AlignmentEntry>();
        addIndel(list, 10, 50, 3, 20); // indel in middle of the read, would be correctly aligned.
        addMutationWithQuality(list, 1, 50, 20, 'A', 'T', Byte.MAX_VALUE);

        Alignments.AlignmentEntry entry = list.get(1);
        assertEquals(20, entry.getSequenceVariations(0).getPosition());
        assertEquals("A", entry.getSequenceVariations(0).getFrom());
        assertEquals("T", entry.getSequenceVariations(0).getTo());
        assertTrue( entry.getSequenceVariations(0).hasToQuality());
        assertEquals(Byte.MAX_VALUE, entry.getSequenceVariations(0).getToQuality().byteAt(0));


    }


    private void addIndel(ObjectList<Alignments.AlignmentEntry> list, int position, int alignLength, int indelLength,
                          int varPosition) {
        Alignments.AlignmentEntry.Builder entry = Alignments.AlignmentEntry.newBuilder();
        MutableString from = new MutableString();
        MutableString to = new MutableString();
        for (int i = 0; i < indelLength; i++) {
            from.append('-');
            to.append('A');
        }
        entry.setTargetIndex(0);
        entry.setMatchingReverseStrand(false);
        entry.setPosition(position);
        entry.setQueryLength(alignLength);
        entry.setQueryAlignedLength(alignLength);
        entry.setQueryIndex(queryIndex++);
        Alignments.SequenceVariation.Builder var = Alignments.SequenceVariation.newBuilder();
        var.setFrom(from.toString());
        var.setTo(to.toString());
        var.setPosition(varPosition);
        var.setReadIndex(varPosition);
        entry.addSequenceVariations(var);
        list.add(entry.build());
    }

    private void addMutation(ObjectList<Alignments.AlignmentEntry> list,
                             int position, int alignLength, int varPosition, char fromBase, char toBase) {
        Alignments.AlignmentEntry.Builder entry = Alignments.AlignmentEntry.newBuilder();
        MutableString from = new MutableString();
        MutableString to = new MutableString();
        from.append(fromBase);
        to.append(toBase);
        entry.setTargetIndex(0);
        entry.setMatchingReverseStrand(false);
        entry.setPosition(position);
        entry.setQueryLength(alignLength);
        entry.setQueryAlignedLength(alignLength);
        entry.setQueryIndex(queryIndex++);
        Alignments.SequenceVariation.Builder var = Alignments.SequenceVariation.newBuilder();
        var.setFrom(from.toString());
        var.setTo(to.toString());
        var.setPosition(varPosition);
        var.setReadIndex(varPosition);
        entry.addSequenceVariations(var);
        list.add(entry.build());
    }

    private void addMutationWithQuality(ObjectList<Alignments.AlignmentEntry> list,
                                        int position, int alignLength, int varPosition, char fromBase, char toBase,
                                        byte toQual) {
        Alignments.AlignmentEntry.Builder entry = Alignments.AlignmentEntry.newBuilder();
        MutableString from = new MutableString();
        MutableString to = new MutableString();
        from.append(fromBase);
        to.append(toBase);
        entry.setTargetIndex(0);
        entry.setMatchingReverseStrand(false);
        entry.setPosition(position);
        entry.setQueryLength(alignLength);
        entry.setQueryAlignedLength(alignLength);
        entry.setQueryIndex(queryIndex++);
        Alignments.SequenceVariation.Builder var = Alignments.SequenceVariation.newBuilder();
        var.setFrom(from.toString());
        var.setTo(to.toString());
        var.setPosition(varPosition);
        var.setReadIndex(varPosition);
        byte[] bytes = new byte[1];
        bytes[0]=toQual;
        final ByteString buffer = ByteString.copyFrom(bytes);
        var.setToQuality(buffer);
        entry.addSequenceVariations(var);
        list.add(entry.build());
    }

}
