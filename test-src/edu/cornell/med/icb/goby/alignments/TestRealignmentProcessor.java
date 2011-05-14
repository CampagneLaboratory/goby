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

import edu.cornell.med.icb.goby.alignments.processors.InfoForTarget;
import edu.cornell.med.icb.goby.alignments.processors.ObservedIndel;
import edu.cornell.med.icb.goby.alignments.processors.RealignmentProcessor;
import edu.cornell.med.icb.goby.modes.AbstractAlignmentToCompactMode;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectListIterator;
import it.unimi.dsi.lang.MutableString;
import static org.junit.Assert.*;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;

/**
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
        assertEquals(2, tinfo.potentialIndels.get(0).getStart());
        assertEquals(5, tinfo.potentialIndels.get(0).getEnd());
        assertEquals(3, tinfo.potentialIndels.get(0).length());


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

                    assertEquals(30 - 1, tinfo.potentialIndels.get(0).getStart());
                    assertEquals(33 - 1, tinfo.potentialIndels.get(0).getEnd());
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
                assertEquals(14, var.getPosition());
                System.out.println("entry:"
                        + entry);
            }

        }
    }

    @Test
    public void testRealignToTheRightReverseStrand() throws IOException {


        ObjectListIterator<Alignments.AlignmentEntry> list = buildList2();
        RealignmentProcessor realigner = new RealignmentProcessor(list);
        realigner.setGenome(new RandomAccessSequenceTestSupport(list2Refs()));
        Alignments.AlignmentEntry entry;
        while ((entry = realigner.nextRealignedEntry(0, 0)) != null) {

            if (entry.getQueryIndex() == 0) {
                Alignments.SequenceVariation var = entry.getSequenceVariations(0);
                assertEquals("CTAG", var.getFrom());
                assertEquals("----", var.getTo());
                assertEquals(14, var.getPosition());
                System.out.println("entry:"
                        + entry);
            }

        }
    }

    @Test
    public void testRealignmentScore1() {
        ObjectListIterator<Alignments.AlignmentEntry> iterator = new ObjectArrayList().iterator();
        RealignmentProcessor realigner = new RealignmentProcessor(iterator);
        RandomAccessSequenceInterface genome = new RandomAccessSequenceTestSupport(list2Refs());
        Alignments.AlignmentEntry entry = makeEntry(0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAATTACTA").build();
        ObservedIndel indel = new ObservedIndel(14, 18, "CTAG", "----");
        assertEquals(3, realigner.score(entry, indel, true, 0, genome));

        entry = makeEntry(0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAATTACTAG").build();
        indel = new ObservedIndel(14, 18, "CTAG", "----");
        assertEquals(4, realigner.score(entry, indel, true, 0, genome));
    }

    @Test
    public void testRealignmentScore2() {
        ObjectListIterator<Alignments.AlignmentEntry> iterator = new ObjectArrayList().iterator();
        RealignmentProcessor realigner = new RealignmentProcessor(iterator);
        RandomAccessSequenceInterface genome = new RandomAccessSequenceTestSupport(list2Refs());
        Alignments.AlignmentEntry entry = makeEntry(0, "ACTGACTGACTGAACTAGTTACTAGCTAAAGTTA", "     CTGACTGAATTAGTA").build();
        ObservedIndel indel = new ObservedIndel(15, 19, "CTAG", "----");
        // expected score is 1 because despite the +3 you get for inserting the indel you have one base mismatch between the ref and the read.
        assertEquals(1, realigner.score(entry, indel, true, 0, genome));
    }

    private int queryIndex = 0;

    private void addEntry(ObjectList<Alignments.AlignmentEntry> list, int targetIndex, String reference, String read) {
        Alignments.AlignmentEntry.Builder entry = makeEntry(targetIndex, reference, read);
        list.add(entry.build());
    }

    private Alignments.AlignmentEntry.Builder makeEntry(int targetIndex, String reference, String read) {
        Alignments.AlignmentEntry.Builder entry = Alignments.AlignmentEntry.newBuilder();

        entry.setTargetIndex(targetIndex);
        entry.setMatchingReverseStrand(false);
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


}
