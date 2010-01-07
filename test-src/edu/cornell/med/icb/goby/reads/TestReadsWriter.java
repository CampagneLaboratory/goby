/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Apr 28, 2009
 *         Time: 12:46:25 PM
 */
public class TestReadsWriter {
    @BeforeClass
    public static void setUp() throws IOException {
        FileUtils.forceMkdir(new File("test-results/reads"));
    }

    @Test
    public void testReadWriteSequences() throws IOException {
        final String[] sequences = {
                "ACTGCGCGCG",
                "AAAAATTTTGGGGGCCCCCCC",
                "AAAAATTTTGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        };
        final String[] descriptions = {
                "hello world",
                "descr2",
                "description 3"
        };

        final ReadsWriter writer = new ReadsWriter(new FileOutputStream(new File("test-results/reads/written-101.bin")));

        int expectedCount = 0;
        writer.setNumEntriesPerChunk(9);
        final int totalNumSeqs = 13;
        for (int i = 0; i < totalNumSeqs; i++) {
            int j = 0;
            for (final String sequence : sequences) {
                writer.setSequence(sequence);
                writer.setDescription(descriptions[j]);
                writer.appendEntry();
                ++j;
                expectedCount++;
            }
        }
        writer.close();
        writer.printStats(System.out);

        final ReadsReader reader = new ReadsReader(new FastBufferedInputStream(new FileInputStream(new File("test-results/reads/written-101.bin"))));
        final MutableString sequence = new MutableString();
        int count = 0;
        while (reader.hasNext()) {
            final Reads.ReadEntry entry = reader.next();
            ReadsReader.decodeSequence(entry, sequence);
            // check that the sequence matches what was encoded:
            assertEquals(sequences[count % sequences.length], sequence.toString());
            assertEquals(descriptions[count % descriptions.length], entry.getDescription());
            assertEquals(count, entry.getReadIndex());
            count++;
        }

        assertEquals(expectedCount, count);
    }

    @Test
    public void testReadFastBufferredOneChunk() throws IOException {
        final String[] sequences = {
                "ACTGCGCGCG",
                "AAAAATTTTGGGGGCCCCCCC",
                "AAAAATTTTGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        };
        final String[] descriptions = {
                "hello world",
                "descr2",
                "description 3"
        };
        final String filename = "test-results/reads/written-101.bin";

        final ReadsWriter writer = new ReadsWriter(new FileOutputStream(new File(filename)));
        ReadsReader reader;

        int expectedCount = 0;
        writer.setNumEntriesPerChunk(9);
        final int totalNumSeqs = 13;
        for (int i = 0; i < totalNumSeqs; i++) {
            int j = 0;
            for (final String sequence : sequences) {
                writer.setSequence(sequence);
                writer.setDescription(descriptions[j]);
                writer.appendEntry();
                ++j;
                expectedCount++;
            }
        }
        writer.close();

        FastBufferedInputStream stream =
                new FastBufferedInputStream(FileUtils.openInputStream(new File(filename)));
        // start at 10 and end at 15 should return no sequences at all.
        reader = new ReadsReader(10, 15, stream);
        MutableString sequence = new MutableString();
        assertFalse(reader.hasNext());

        stream = new FastBufferedInputStream(FileUtils.openInputStream(new File(filename)));
        // start at 138 (precise start of a chunk) and end at 139 should return 9 sequences exactly.
        reader = new ReadsReader(138, 139, stream);
        sequence = new MutableString();
        int count = 0;

        while (reader.hasNext()) {
            final Reads.ReadEntry entry = reader.next();
            ReadsReader.decodeSequence(entry, sequence);
            // check that the sequence matches what was encoded:
            assertEquals(sequences[count % sequences.length], sequence.toString());
            assertEquals(descriptions[count % descriptions.length], entry.getDescription());
            //  assertEquals(count+9, entry.getReadIndex());
            //        System.out.println("desc: " +entry.getDescription() );
            count++;
        }
        // should have skipped 2:
        assertEquals(9, count);
    }

    @Test
    public void testReadFastBufferedInChunks() throws IOException {
        final String[] sequences = {
                "ACTGCGCGCG",
                "AAAAATTTTGGGGGCCCCCCC",
                "AAAAATTTTGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        };
        final String[] descriptions = {
                "hello world",
                "descr2",
                "description 3"
        };
        final String filename = "test-results/reads/written-101.bin";
        final ReadsWriter writer = new ReadsWriter(new FileOutputStream(new File(filename)));

        int expectedCount = 0;
        writer.setNumEntriesPerChunk(9);
        final int totalNumSeqs = 13;
        for (int i = 0; i < totalNumSeqs; i++) {
            int j = 0;
            for (final String sequence : sequences) {
                writer.setSequence(sequence);
                writer.setDescription(descriptions[j]);
                writer.appendEntry();
                ++j;
                expectedCount++;
            }
        }
        writer.close();

        FastBufferedInputStream stream =
                new FastBufferedInputStream(FileUtils.openInputStream(new File(filename)));
        writer.printStats(System.out);

        final MutableString sequence = new MutableString();
        int count = 0;

        ReadsReader reader = new ReadsReader(0, Long.MAX_VALUE, stream);

        while (reader.hasNext()) {
            final Reads.ReadEntry entry = reader.next();
            ReadsReader.decodeSequence(entry, sequence);
            // check that the sequence matches what was encoded:
            assertEquals(sequences[count % sequences.length], sequence.toString());
            assertEquals(descriptions[count % descriptions.length], entry.getDescription());
            assertEquals(count, entry.getReadIndex());
            count++;
        }

        assertEquals(expectedCount, count);

        stream.close();

        // now start at position 10, will skip the first collection.
        stream = new FastBufferedInputStream(FileUtils.openInputStream(new File(filename)));
        reader = new ReadsReader(10, Long.MAX_VALUE, stream);

        count = 0;
        while (reader.hasNext()) {
            final Reads.ReadEntry entry = reader.next();
            ReadsReader.decodeSequence(entry, sequence);
            // check that the sequence matches what was encoded:
            assertEquals(sequences[count % sequences.length], sequence.toString());
            assertEquals(descriptions[count % descriptions.length], entry.getDescription());
            //  assertEquals(count+9, entry.getReadIndex());
            count++;
        }

        assertEquals(expectedCount - 9, count);

        // start at 150 is past the start of the second read collection, so this collection is not returned:
        reader = new ReadsReader(150, Long.MAX_VALUE, stream);

        count = 0;
        while (reader.hasNext()) {
            final Reads.ReadEntry entry = reader.next();
            ReadsReader.decodeSequence(entry, sequence);
            // check that the sequence matches what was encoded:
            assertEquals(sequences[count % sequences.length], sequence.toString());
            assertEquals(descriptions[count % descriptions.length], entry.getDescription());
            //  assertEquals(count+9, entry.getReadIndex());
            count++;
        }
        // should have skipped 2:
        assertEquals(expectedCount - 9 * 2, count);

        // start at 0 and end at 150 should return only the second read collection, 9 sequences exactly.
        reader = new ReadsReader(10, 150, stream);

        count = 0;
        while (reader.hasNext()) {
            final Reads.ReadEntry entry = reader.next();
            ReadsReader.decodeSequence(entry, sequence);
            // check that the sequence matches what was encoded:
            assertEquals(sequences[count % sequences.length], sequence.toString());
            assertEquals(descriptions[count % descriptions.length], entry.getDescription());
            //  assertEquals(count+9, entry.getReadIndex());
            count++;
        }

        // should have skipped 2:
        assertEquals(9, count);
    }
}
