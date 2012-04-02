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

import edu.cornell.med.icb.goby.readers.FastXEntry;
import edu.cornell.med.icb.goby.readers.FastXReader;
import edu.cornell.med.icb.goby.reads.QualityEncoding;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceCache;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.bytes.ByteList;
import it.unimi.dsi.fastutil.chars.CharArrayList;
import it.unimi.dsi.fastutil.chars.CharList;
import it.unimi.dsi.lang.MutableString;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * Test ExportableAlignmentEntryData, a class which assists with exporting from Compact to other formats such
 * as SAM.
 */
public class TestExportableAlignmentEntryData {

    @Test
    public void testRandomAccess() throws IOException, ClassNotFoundException {
        final RandomAccessSequenceCache genome = new RandomAccessSequenceCache();
        genome.load("test-data/seq-var-test/small-synth.random-access-genome");
        
        final List<Character> refBases = new ArrayList<Character>();
        for (int i = 10; i <= 15; i++) {
            final char base = genome.get(0, i);
            refBases.add(base);
        }
        assertTrue("Incorrect base read from randomSequenceCache", 'G' == refBases.get(0));
        assertTrue("Incorrect base read from randomSequenceCache", 'C' == refBases.get(1));
        assertTrue("Incorrect base read from randomSequenceCache", 'T' == refBases.get(2));
        assertTrue("Incorrect base read from randomSequenceCache", 'G' == refBases.get(3));
        assertTrue("Incorrect base read from randomSequenceCache", 'G' == refBases.get(4));
        assertTrue("Incorrect base read from randomSequenceCache", 'A' == refBases.get(5));
    }

    @Test
    public void testMakeExportable() throws IOException, ClassNotFoundException {
        final FastXReader fastqReader = new FastXReader("test-data/seq-var-test/seq-var-reads.fq");
        final Map<Integer, ReadsDataEntry> reads = new LinkedHashMap<Integer, ReadsDataEntry>();
        int readIndex = 0;
        for (final FastXEntry fastqEntry : fastqReader) {
            final ReadsDataEntry entry = new ReadsDataEntry(readIndex++,
                    fastqEntry.getEntryHeader(),
                    fastqEntry.getSequence(),
                    fastqEntry.getQuality(), -64);
            reads.put(entry.readIndex, entry);
        }


        final AlignmentReader reader = new AlignmentReaderImpl("test-data/seq-var-test/sorted-seq-var-reads-gsnap.entries");
        final RandomAccessSequenceCache genome = new RandomAccessSequenceCache();
        genome.load("test-data/seq-var-test/small-synth.random-access-genome");
        final ExportableAlignmentEntryData exportData = new ExportableAlignmentEntryData(genome, QualityEncoding.PHRED);
        while (reader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.next();
            final ReadsDataEntry actualReadsEntry = reads.get(alignmentEntry.getQueryIndex());
            System.out.printf("Processing queryIndex=%d with description '%s'%n", alignmentEntry.getQueryIndex(),
                    actualReadsEntry.readName);
            exportData.buildFrom(alignmentEntry, actualReadsEntry.readBases, actualReadsEntry.readQuals);
            assertFalse(exportData.getInvalidMessage(), exportData.isInvalid());
        }
    }

    private class ReadsDataEntry {
        final int readIndex;
        final String readName;
        final CharList readBases;
        final ByteList readQuals;

        private ReadsDataEntry(final int readIndex, final MutableString readName,
                              final MutableString readBases, final MutableString readQuals, final int qualityOffset) {
            this.readIndex = readIndex;
            this.readName = readName.toString();
            this.readBases = new CharArrayList(readBases.length());
            for (int i = 0; i < readBases.length(); i++) {
                this.readBases.add(readBases.charAt(i));
            }
            this.readQuals = new ByteArrayList(readQuals.length());
            for (int i = 0; i < readQuals.length(); i++) {
                this.readQuals.add((byte) ((byte) readQuals.charAt(i) + qualityOffset));
            }
        }
    }
}
