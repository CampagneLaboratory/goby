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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.lang.MutableString;
import org.junit.Test;

import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;

import static junit.framework.Assert.assertEquals;


/**
 * @author Fabien Campagne
 *         Date: 9/14/12
 *         Time: 11:05 AM
 */
public class ExportSplicingEventsTest {
    @Test
    public void testProcess() throws Exception {
        final Collection<Alignments.AlignmentEntry> collection = new ArrayList<Alignments.AlignmentEntry>();
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 34)).build());
        String expected = "sample-id\t1\t31\t35\t+\t??\t1\t0.0\n";

        assertConditions(collection, expected);
    }

    @Test
    public void testProcess2() throws Exception {
        final Collection<Alignments.AlignmentEntry> collection = new ArrayList<Alignments.AlignmentEntry>();
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 34)).build());
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 34)).build());
        String expected = "sample-id\t1\t31\t35\t+\t??\t2\t0.0\n";

        assertConditions(collection, expected);
    }

    @Test
    public void testProcess3() throws Exception {
        final Collection<Alignments.AlignmentEntry> collection = new ArrayList<Alignments.AlignmentEntry>();
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setTargetAlignedLength(3)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 34)).build());
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 33)).build());         // end off by one base
        String expected = "sample-id\t1\t31\t34\t+\t??\t1\t0.0\n" +
                "sample-id\t1\t34\t35\t+\t??\t1\t0.0\n";
        assertConditions(collection, expected);
    }


    @Test
    public void testProcess4() throws Exception {
        final Collection<Alignments.AlignmentEntry> collection = new ArrayList<Alignments.AlignmentEntry>();
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 34)).build());
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 33)).build());         // end off by one base
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 34)).build());
        String expected = "sample-id\t1\t31\t35\t+\t??\t2\t0.0\n" +
                "sample-id\t1\t31\t34\t+\t??\t1\t0.0\n";

        assertConditions(collection, expected);
    }

    @Test
    public void testProcessSort() throws Exception {
        final Collection<Alignments.AlignmentEntry> collection = new ArrayList<Alignments.AlignmentEntry>();
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setTargetAlignedLength(10)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 34)).build());
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setTargetAlignedLength(9)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 33)).build());         // end off by one base
        collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                .setQueryLength(40)
                .setPosition(30)
                .setTargetAlignedLength(5)
                .setMatchingReverseStrand(false)
                .setSplicedForwardAlignmentLink(newSpliceLink(0, 34)).build());
        String expected = "sample-id\t1\t36\t35\t+\t??\t1\t0.0\n" +
                "sample-id\t1\t40\t34\t+\t??\t1\t0.0\n" +
                "sample-id\t1\t41\t35\t+\t??\t1\t0.0\n";

        assertConditions(collection, expected);
    }
    @Test
       public void testMappingQualityThreshold() throws Exception {
           final Collection<Alignments.AlignmentEntry> collection = new ArrayList<Alignments.AlignmentEntry>();
           collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                   .setQueryLength(40)
                   .setPosition(30)
                   .setMappingQuality(255)
                   .setMatchingReverseStrand(false)
                   .setSplicedForwardAlignmentLink(newSpliceLink(0, 34)).build());
           collection.add(Alignments.AlignmentEntry.newBuilder().setTargetIndex(0)
                   .setQueryLength(40)
                   .setPosition(30)
                   .setMappingQuality(254)
                   .setMatchingReverseStrand(false)
                   .setSplicedForwardAlignmentLink(newSpliceLink(0, 34)).build());
           String expected = "sample-id\t1\t31\t35\t+\t??\t1\t0.0\n";

           assertConditions(collection, expected);
       }
    private void assertConditions(Collection<Alignments.AlignmentEntry> collection, String expected) throws IOException {
        final StringWriter output = new StringWriter();
        final ExportSplicingEvents processor = new ExportSplicingEvents(output);
        final IndexedIdentifier ids = new IndexedIdentifier();
        ids.registerIdentifier(new MutableString("1"));
        ids.registerIdentifier(new MutableString("2"));
        final DoubleIndexedIdentifier reverseIds = new DoubleIndexedIdentifier(ids);

        processor.process(collection, reverseIds);
        assertEquals(expected, output.getBuffer().toString());
    }

    private Alignments.RelatedAlignmentEntry newSpliceLink(final int targetIndex, final int end) {
        return Alignments.RelatedAlignmentEntry.newBuilder().setTargetIndex(targetIndex).setPosition(end).build();
    }
}
