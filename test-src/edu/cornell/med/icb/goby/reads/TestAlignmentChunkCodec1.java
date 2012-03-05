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
        final AlignmentReader reader = new AlignmentReaderImpl("/data/rrbs/AJPBRWE.entries");
        final Alignments.AlignmentCollection.Builder collectionBuilder = Alignments.AlignmentCollection.newBuilder();
       int i=0;
        for (Alignments.AlignmentEntry entry : reader) {
            collectionBuilder.addAlignmentEntries(entry);
            if (i++>10000) {
                break;
            }
        }

        final ByteArrayOutputStream encoded = codec.encode(collectionBuilder.build());
        assertNotNull(encoded);
       // assertEquals(collectionBuilder.getAlignmentEntriesCount(), encoded.);
    }
}
