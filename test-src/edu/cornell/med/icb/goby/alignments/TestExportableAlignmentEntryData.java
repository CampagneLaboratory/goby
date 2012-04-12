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
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

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
        final RandomAccessSequenceCache genome = new RandomAccessSequenceCache();
        genome.load("test-data/seq-var-test/small-synth.random-access-genome");
        Int2ObjectMap<Map<String, String>> samDetailsMap = readSamFileToMap(genome, "test-data/seq-var-test/new/seq-var-reads-gsnap.sam");

        final FastXReader fastqReader = new FastXReader("test-data/seq-var-test/new/seq-var-reads.fq");
        final Map<Integer, ReadsDataEntry> reads = new LinkedHashMap<Integer, ReadsDataEntry>();
        int readIndex = 0;
        for (final FastXEntry fastqEntry : fastqReader) {
            final ReadsDataEntry entry = new ReadsDataEntry(readIndex++,
                    fastqEntry.getEntryHeader(),
                    fastqEntry.getSequence(),
                    fastqEntry.getQuality(), -64);
            reads.put(entry.readIndex, entry);
        }

        final AlignmentReader reader = new AlignmentReaderImpl("test-data/seq-var-test/new/sorted-seq-var-reads-gsnap.entries");
        final ExportableAlignmentEntryData exportData = new ExportableAlignmentEntryData(genome, QualityEncoding.PHRED);
        while (reader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.next();
            final ReadsDataEntry actualReadsEntry = reads.get(alignmentEntry.getQueryIndex());
            if (actualReadsEntry == null) {
                fail("Couldn't find actual read for alignmentEntry.queryIndex=" + alignmentEntry.getQueryIndex());
            }
            System.out.printf("Processing queryIndex=%d with description '%s'%n", alignmentEntry.getQueryIndex(),
                    actualReadsEntry.readName);
            exportData.buildFrom(alignmentEntry, actualReadsEntry.readBases, actualReadsEntry.readQuals);
            assertFalse(exportData.getInvalidMessage(), exportData.isInvalid());
            validateEntry(exportData, samDetailsMap.get(exportData.getQueryIndex()));
        }
    }

    private void validateEntry(final ExportableAlignmentEntryData exportData,
                                  final Map<String, String> samData) {
        final int qi = exportData.getQueryIndex();
        if (samData == null) {
            fail("Couldn't find sam data for queryIndex=" + qi);
        }
        assertEquals("Pair flags are incorrect for qi=" + qi, (int) Integer.valueOf(samData.get("pairFlags")), exportData.getPairFlags());
        assertEquals("targetIndex incorrect for qi=" + qi, (int) Integer.valueOf(samData.get("targetIndex")), exportData.getTargetIndex());
        assertEquals("position incorrect for qi=" + qi, (int) Integer.valueOf(samData.get("position")), exportData.getStartPosition());
        assertEquals("mapq incorrect for qi=" + qi, (int) Integer.valueOf(samData.get("mapq")), exportData.getMappingQuality());
        assertEquals("cigar incorrect for qi=" + qi, samData.get("cigar"), exportData.getCigarString());
        validateSequence(qi, samData.get("read"), exportData.getReadBasesOriginal());
        assertEquals("mismatches incorrect for qi=" + qi, samData.get("mismatches"), "MD:Z:" + exportData.getMismatchString());
    }

    private void validateSequence(final int qi, final String expected, final String actual) {
        assertEquals("Read lengths don't match qi=" + qi, expected.length(), actual.length());
        for (int i = 0; i < expected.length(); i++) {
            final char actualBase = actual.charAt(i);
            final char expectedBase = expected.charAt(i);
            if (actualBase != expectedBase && actualBase != 'N') {
                assertEquals("Base incorrect at qi=" + qi + " i=" + i, expectedBase, actualBase);
            }
        }
    }


    private Int2ObjectMap<Map<String, String>> readSamFileToMap(final RandomAccessSequenceCache genome, final String filename) throws IOException {
        List<String> lines = FileUtils.readLines(new File(filename));
        Int2ObjectMap<Map<String, String>> result = new Int2ObjectArrayMap<Map<String, String>>();
        for (final String line : lines) {
            if (line.startsWith("@")) {
                continue;
            }
            final String[] parts = line.split("\t");
            final String target = parts[2];
            if (target.equals("*")) {
                continue;
            }
            final int queryIndex = Integer.valueOf(parts[0]);
            Map<String, String> entry = new HashMap<String, String>();
            result.put(queryIndex, entry);
            entry.put("pairFlags", parts[1]);
            entry.put("targetIndex", String.valueOf(genome.getReferenceIndex(target)));
            entry.put("position", parts[3]);
            entry.put("mapq", parts[4]);
            entry.put("cigar", parts[5]);
            entry.put("read", parts[9]);
            for (int i = 10; i < parts.length; i++) {
                if (parts[i].startsWith("MD:Z")) {
                    entry.put("mismatches", parts[i]);
                    break;
                }
            }
        }
        return result;
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
