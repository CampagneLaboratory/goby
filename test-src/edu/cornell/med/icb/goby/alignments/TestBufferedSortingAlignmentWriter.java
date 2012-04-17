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

package edu.cornell.med.icb.goby.alignments;

import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 4/17/12
 *         Time: 5:47 PM
 */
public class TestBufferedSortingAlignmentWriter {

    private String expectedCase1 = "{target_index: 0\n" +
            "position: 12\n" +
            "}\n" +
            "{target_index: 0\n" +
            "position: 13\n" +
            "}\n" +
            "{target_index: 0\n" +
            "position: 14\n" +
            "}\n" +
            "{target_index: 0\n" +
            "position: 15\n" +
            "}\n" +
            "Closed\n";



    @Test
    // ordered entries
    public void testCase1() throws IOException {
        AlignmentToTextWriter destination = new AlignmentToTextWriter();
        BufferedSortingAlignmentWriter writer = new BufferedSortingAlignmentWriter(destination);
        final Alignments.AlignmentEntry entry1 = buildEntryWithTargetPosition(0, 12);
        final Alignments.AlignmentEntry entry2 = buildEntryWithTargetPosition(0, 13);
        final Alignments.AlignmentEntry entry3 = buildEntryWithTargetPosition(0, 14);
        final Alignments.AlignmentEntry entry4 = buildEntryWithTargetPosition(0, 15);
        writer.appendEntry(entry1);
        writer.appendEntry(entry2);
        writer.appendEntry(entry3);
        writer.appendEntry(entry4);
        writer.close();
        assertEquals(expectedCase1, destination.getTextOutput().toString());
    }

    @Test
    // ordered entries
    public void testCase2() throws IOException {
        AlignmentToTextWriter destination = new AlignmentToTextWriter();
        BufferedSortingAlignmentWriter writer = new BufferedSortingAlignmentWriter(destination);
        final Alignments.AlignmentEntry entry1 = buildEntryWithTargetPosition(0, 12);
        final Alignments.AlignmentEntry entry2 = buildEntryWithTargetPosition(0, 13);
        final Alignments.AlignmentEntry entry3 = buildEntryWithTargetPosition(0, 14);
        final Alignments.AlignmentEntry entry4 = buildEntryWithTargetPosition(0, 15);
        writer.appendEntry(entry4);
        writer.appendEntry(entry2);
        writer.appendEntry(entry3);
        writer.appendEntry(entry1);
        writer.close();
        assertEquals(expectedCase1, destination.getTextOutput().toString());
    }

    @Test
    // ordered entries
    public void testCase3() throws IOException {
        AlignmentToTextWriter destination = new AlignmentToTextWriter();
        BufferedSortingAlignmentWriter writer = new BufferedSortingAlignmentWriter(destination);
        final Alignments.AlignmentEntry entry1 = buildEntryWithTargetPosition(5, 12);
        final Alignments.AlignmentEntry entry2 = buildEntryWithTargetPosition(4, 13);
        final Alignments.AlignmentEntry entry3 = buildEntryWithTargetPosition(3, 14);
        final Alignments.AlignmentEntry entry4 = buildEntryWithTargetPosition(2, 15);
        writer.appendEntry(entry4);
        writer.appendEntry(entry2);
        writer.appendEntry(entry3);
        writer.appendEntry(entry1);
        writer.close();
        assertEquals(expectedCase3, destination.getTextOutput().toString());
    }
    private String expectedCase3 = "{target_index: 2\n" +
            "position: 15\n" +
            "}\n" +
            "{target_index: 3\n" +
            "position: 14\n" +
            "}\n" +
            "{target_index: 4\n" +
            "position: 13\n" +
            "}\n" +
            "{target_index: 5\n" +
            "position: 12\n" +
            "}\n" +
            "Closed\n";


    private Alignments.AlignmentEntry buildEntryWithTargetPosition(int targetIndex, int position) {
        Alignments.AlignmentEntry.Builder builder = Alignments.AlignmentEntry.newBuilder();
        builder.setTargetIndex(targetIndex).setPosition(position);
        return builder.build();
    }
}
