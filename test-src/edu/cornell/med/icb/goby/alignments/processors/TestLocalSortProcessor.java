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

package edu.cornell.med.icb.goby.alignments.processors;

import edu.cornell.med.icb.goby.alignments.Alignments;
import org.easymock.EasyMock;
import org.junit.Before;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNull;
import static org.easymock.EasyMock.*;

/**
 * @author Fabien Campagne
 *         Date: 6/2/11
 *         Time: 12:00 PM
 */
public class TestLocalSortProcessor {
    AlignmentProcessorInterface delegate;
    private LocalSortProcessor resortProcessor;
    private Alignments.AlignmentEntry entry;

    @Before
    public void setUp() throws Exception {
        delegate = EasyMock.createMock(AlignmentProcessorInterface.class);
        resortProcessor = new LocalSortProcessor(delegate);
    }

    @Test
    public void testEmpty() throws Exception {
        replay();
        assertNull(resortProcessor.nextRealignedEntry(0, 0));
    }


    @Test
    public void testOneElement() throws Exception {
        entry = buildEntry(1, 0, 1);
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry);
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(null);
        replay(delegate);
        assertEquals(entry, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals(null, resortProcessor.nextRealignedEntry(0, 0));
        verify(delegate);
    }

    @Test
    public void testTwoElements() throws Exception {

        Alignments.AlignmentEntry entry_0_10;
        Alignments.AlignmentEntry entry_0_1;
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_0_10 = buildEntry(2, 0, 10));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_0_1 = buildEntry(1, 0, 1));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(null);
        replay(delegate);
        assertEquals("entry_0_1 must have been reordered first", entry_0_1, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("entry_0_10 must have been reordered second", entry_0_10, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("", null, resortProcessor.nextRealignedEntry(0, 0));
        verify(delegate);
    }

    @Test
    public void testSeveralElements() throws Exception {

        Alignments.AlignmentEntry entry_0_9;
        Alignments.AlignmentEntry entry_0_10;
        Alignments.AlignmentEntry entry_0_14;
        Alignments.AlignmentEntry entry_0_13;
        Alignments.AlignmentEntry entry_0_12;
        Alignments.AlignmentEntry entry_0_11;
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_0_9 = buildEntry(2, 0, 9));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_0_9 = buildEntry(2, 0, 9));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_0_11 = buildEntry(2, 0, 11));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_0_12 = buildEntry(2, 0, 12));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_0_10 = buildEntry(1, 0, 10));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_0_13 = buildEntry(2, 0, 13));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_0_14 = buildEntry(2, 0, 14));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(null);
        replay(delegate);
        assertEquals("entry_0_9 must have been reordered first", entry_0_9, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("entry_0_9 must have been reordered second", entry_0_9, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("entry_0_10 must have been reordered third", entry_0_10, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("entry_0_11 must have been reordered fourth", entry_0_11, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("entry_0_12 must have been reordered fith", entry_0_12, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("entry_0_13 must have been reordered sixth", entry_0_13, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("entry_0_14 must have been reordered seventh", entry_0_14, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("", null, resortProcessor.nextRealignedEntry(0, 0));
        verify(delegate);
    }

    @Test
    public void testTwoTargets() throws Exception {

        Alignments.AlignmentEntry entry_0_10;
        Alignments.AlignmentEntry entry_1_1;
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_0_10 = buildEntry(1, 0, 10));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(entry_1_1 = buildEntry(2, 1, 1));
        expect(delegate.nextRealignedEntry(0, 0)).andReturn(null);
        replay(delegate);
        assertEquals("expected entry on target 0", entry_0_10, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("expected entry on target 1", entry_1_1, resortProcessor.nextRealignedEntry(0, 0));
        assertEquals("", null, resortProcessor.nextRealignedEntry(0, 0));
        verify(delegate);
    }

    private Alignments.AlignmentEntry buildEntry(int queryIndex, int targetIndex, int position) {
        return Alignments.AlignmentEntry.newBuilder().
                setQueryIndex(queryIndex).
                setTargetIndex(targetIndex).
                setPosition(position).
                setMatchingReverseStrand(false).build();
    }
}
