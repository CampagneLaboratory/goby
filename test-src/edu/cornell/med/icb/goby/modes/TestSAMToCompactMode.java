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

package edu.cornell.med.icb.goby.modes;

import org.junit.Test;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import edu.cornell.med.icb.goby.alignments.PerQueryAlignmentData;
import edu.cornell.med.icb.goby.alignments.TestIteratedSortedAlignment2;
import edu.cornell.med.icb.goby.alignments.IterateSortedAlignmentsTester;

/**
 * Created by IntelliJ IDEA.
 * User: kdorff
 * Date: Apr 13, 2011
 * Time: 12:27:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class TestSAMToCompactMode {

    Int2ObjectMap<PerQueryAlignmentData> alignmentDataMap;
    int[] alignmentQueryIndexes;
    Int2ObjectMap<PerQueryAlignmentData> seqvarDataMap;
    int[] seqvarQueryIndexes;
    boolean dataSetup = false;

    /**
     * Test that the queryIndexes from the alignment (ONLY those queryIndexes that have sequence variations)
     * match with those from the "seq-var-reads-gsnap.display-seq-var-tsv-base.tsv" file.
     * @throws IOException error reading input files
     */
    @Test
    public void testIndexesMatch() throws IOException {
        setupData();
        assertArrayEquals("alignmentDataMap indexes doesn't match seqvarDataMap indexes",
                seqvarQueryIndexes, alignmentQueryIndexes);
    }

    /**
     * Test that the sequence variations in the alignment match those we know to be good from
     * the "seq-var-reads-gsnap.display-seq-var-tsv-base.tsv" file.
     * @throws IOException error reading input files
     */
    @Test
    public void testSequenceVariationsMatch() throws IOException {
        setupData();
        for (int queryIndex : seqvarQueryIndexes) {
            PerQueryAlignmentData align = alignmentDataMap.get(queryIndex);
            PerQueryAlignmentData var = seqvarDataMap.get(queryIndex);

            Map<String, String> alignSeqVarsMap = align.refPositionReadIndexToBaseMap;
            Map<String, String> varSeqVarsMap = var.refPositionReadIndexToBaseMap;

            assertEquals(String.format("queryIndex=%d alignSeqVarsMap.size()(%d) should equal varSeqVarsMap.size()(%d)",
                    queryIndex,
                    alignSeqVarsMap.size(), varSeqVarsMap.size()),
                    alignSeqVarsMap.size(), varSeqVarsMap.size());
            for (Map.Entry<String, String> varEntry : varSeqVarsMap.entrySet()) {
                // Make sure the sequence variations match
                final String varEntryBases = varEntry.getValue();
                final String alignEntryBases = alignSeqVarsMap.get(varEntry.getKey());
                assertNotNull(String.format("queryIndex=%d Could not find alignSeqVarsMap entry for %s",
                        queryIndex, varEntry.getKey()),
                        alignEntryBases);
                assertEquals(String.format("queryIndex=%d alignEntryBases(%s) should equal varEntryBases(%s)",
                        queryIndex, alignEntryBases, varEntryBases),
                        alignEntryBases, varEntryBases);
            }
        }
    }

    public synchronized void setupData() throws IOException {
        if (dataSetup) {
            return;
        }
        System.out.println("++");
        System.out.println("++ Running setupData");
        System.out.println("++");
        SAMToCompactOldMode converter = new SAMToCompactOldMode();
        converter.setInputFile("test-data/sam-to-compact-test/seq-var-reads-gsnap.sam");
        final String testOutputDir = "test-results/sam-to-compact-test";
        final String testOutputBasename = "seq-var-reads-gsnap";
        final File testOutputDirFile = new File(testOutputDir);
        FileUtils.deleteDirectory(testOutputDirFile);
        testOutputDirFile.mkdirs();
        converter.setOutputFile(testOutputDir + "/presort-" + testOutputBasename);
        converter.setQueryReadIdsFilename("test-data/sam-to-compact-test/seq-var-reads.compact-reads");
        // For this test, we want to keep ALL records
        converter.setQualityFilterParameters("threshold=1");
        // converter.setPropagateQueryIds(true);
        converter.setPropagateTargetIds(true);
        converter.execute();

        // Sort the data
        SortMode sorter = new SortMode();
        sorter.setInput(testOutputDir + "/presort-" + testOutputBasename);
        sorter.setOutput(testOutputDir + "/" + testOutputBasename);
        sorter.execute();

        // Read the seq-vars from the new sam file
        IterateSortedAlignmentsTester alignmentIterator = new IterateSortedAlignmentsTester();
        final String[] singleBasename = {testOutputDir + "/" + testOutputBasename};
        alignmentIterator.iterate(singleBasename);
        alignmentIterator.removeWithoutSeqvars();
        alignmentDataMap = alignmentIterator.queryIndexToAlignmentDataMap;
        alignmentQueryIndexes = alignmentDataMap.keySet().toIntArray();
        Arrays.sort(alignmentQueryIndexes);

        // Read the seq-vars from the known good file
        seqvarDataMap = TestIteratedSortedAlignment2.readSeqVarFile(
                "test-data/seq-var-test/seq-var-reads-gsnap.display-seq-var-tsv-base.tsv");
        seqvarQueryIndexes = seqvarDataMap.keySet().toIntArray();
        Arrays.sort(seqvarQueryIndexes);
        dataSetup = true;
    }
}
