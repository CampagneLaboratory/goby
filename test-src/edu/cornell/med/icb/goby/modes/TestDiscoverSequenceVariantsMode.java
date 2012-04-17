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

import com.google.protobuf.ByteString;
import com.martiansoftware.jsap.JSAPException;
import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.util.TestFiles;
import edu.cornell.med.icb.io.TSVReader;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;
import java.util.Collections;

import static org.junit.Assert.assertNotNull;

/**
 * @author Fabien Campagne
 *         Date: Mar 12, 2011
 *         Time: 12:15:11 PM
 */
public class TestDiscoverSequenceVariantsMode extends TestFiles {


    private static final Log LOG = LogFactory.getLog(TestDiscoverSequenceVariantsMode.class);
    private static final String BASE_TEST_DIR = "test-results/discover-variants";
    private static final int NUM_BASENAME_PER_GROUP = 5;
    private static final String[] basenames = new String[NUM_BASENAME_PER_GROUP * 2];
    private static final int NUM_TRIALS = 10;
    private static String statsFilename = BASE_TEST_DIR + "/variations-stats2.tsv";
    private static String[] specificAlignments;


    @BeforeClass
    public static void init() {

        GobyRengine.getInstance();
        //    assertNotNull("R engine must be available", GobyRengine.getInstance().getRengine());
    }

    @Test
    public void testDiscoverGroupsOnly() throws IOException, JSAPException {

        DiscoverSequenceVariantsMode mode = new DiscoverSequenceVariantsMode();
        // loop to check that results are reproducible and do not change from execution to execution..
        for (int i = 0; i < NUM_TRIALS; i++) {
            String outputFilename = "out" + i + ".tsv";
            String[] args = constructArgumentString(
                    basenames, BASE_TEST_DIR + "/" + outputFilename, "none").split("[\\s]");
            args = add(args, new String[]{"--format", DiscoverSequenceVariantsMode.OutputFormat.BETWEEN_GROUPS.toString()});


            configureTestGenome(mode);
            mode.configure(args);
            mode.setCallIndels(false);
            mode.setDisableAtLeastQuarterFilter(true);
            mode.execute();

            assertEquals(new File(BASE_TEST_DIR + "/" + outputFilename),
                    new File("test-data/discover-variants/expected-groups-only.tsv"));

        }
    }

    private void configureTestGenome(DiscoverSequenceVariantsMode mode) {
        RandomAccessSequenceTestSupport genomeForTest = new RandomAccessSequenceTestSupport(new String[]{"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"}) {
            @Override
            public char get(int referenceIndex, int position) {

                return 'C';
            }
        };
        mode.setTestGenome(genomeForTest);
    }

    @Test
    public void testDiscoverWithSamples() throws IOException, JSAPException {

        DiscoverSequenceVariantsMode mode = new DiscoverSequenceVariantsMode();
        int i = 1;
        String outputFilename = "out-samples-" + i + ".tsv";
        String[] args = constructArgumentString(
                basenames, BASE_TEST_DIR + "/" + outputFilename, "samples").split("[\\s]");
        args = add(args, new String[]{"--format", DiscoverSequenceVariantsMode.OutputFormat.BETWEEN_GROUPS.toString()});

        configureTestGenome(mode);
        mode.configure(args);
        mode.setCallIndels(false);
        mode.setDisableAtLeastQuarterFilter(true);
        mode.execute();
        assertEquals(new File(BASE_TEST_DIR + "/" + outputFilename),
                new File("test-data/discover-variants/expected-output-samples.tsv"));

    }


    //not a test, just useful to run through the code.

    public void testDiscoverMethylationOneGroup() throws IOException, JSAPException {

        final DiscoverSequenceVariantsMode mode = new DiscoverSequenceVariantsMode();
        int i = 1;
        final String outputFilename = "out-methylation-" + i + ".tsv";
        String[] args = constructArgumentString(
                basenames, BASE_TEST_DIR + "/" + outputFilename, "methylation").split("[\\s]");
        args = add(args, new String[]{"--format", DiscoverSequenceVariantsMode.OutputFormat.METHYLATION.toString()});
        for (int j = 0; j < args.length; j++) {
            if ("--compare".equals(args[j])) {
                // remove the compare argument:
                args[j + 1] = null;
                args[j] = null;
            }
        }
        args = removeNulled(args);
        configureTestGenome(mode);
        mode.configure(args);
        mode.execute();
        assertEquals(new File(BASE_TEST_DIR + "/" + outputFilename),
                new File("test-data/discover-variants/expected-output-methylation.tsv"));

    }

    private String[] removeNulled(final String[] args) {
        int removeCount = 0;
        for (final String s : args) {
            if (s == null) {
                removeCount++;
            }
        }
        final String[] dest = new String[args.length - removeCount];
        int j = 0;
        for (final String s : args) {
            if (s != null) {
                dest[j++] = s;
            }
        }
        return dest;
    }

    @Test
    public void testDiscoverAlleleImbalance() throws IOException, JSAPException {

        DiscoverSequenceVariantsMode mode = new DiscoverSequenceVariantsMode();
        int i = 1;
        String outputFilename = "out-allele-" + i + ".tsv";
        String[] args = constructArgumentString(
                basenames, BASE_TEST_DIR + "/" + outputFilename, "samples").split("[\\s]");
        args = add(args, new String[]{"--format", DiscoverSequenceVariantsMode.OutputFormat.ALLELE_FREQUENCIES.toString()});

        configureTestGenome(mode);
        FormatConfigurator configurator = new FormatConfigurator() {

            @Override
            public void configureFormatter(final SequenceVariationOutputFormat formatter) {
                ((AlleleFrequencyOutputFormat) formatter).setMinimumAllelicDifference(0);
            }
        };
        mode.setFormatConfigurator(configurator);
        mode.configure(args);

        mode.execute();
        assertEquals(
                new File("test-data/discover-variants/expected-output-alleles.tsv"),
                new File(BASE_TEST_DIR + "/" + outputFilename)
        );

    }

    @Test
    public void testDiscoverCheckSampleCountAssociation() throws IOException, JSAPException {

        DiscoverSequenceVariantsMode mode = new DiscoverSequenceVariantsMode();
        // loop to check that results are reproducible and do not change from execution to execution..
        for (int i = 0; i < NUM_TRIALS; i++) {
            String outputFilename = "out-sample-counts" + i + ".tsv";
            final String[] basenameArray = shuffle(add(basenames, specificAlignments));
            final String output = BASE_TEST_DIR + "/" + outputFilename;
            String[] args = constructArgumentString(
                    basenameArray, output, "samples").split("[\\s]");
            args = add(args, new String[]{"--format", DiscoverSequenceVariantsMode.OutputFormat.BETWEEN_GROUPS.toString()});

            configureTestGenome(mode);
            mode.configure(args);

            mode.execute();

            TSVReader reader = new TSVReader(new FileReader(output));
            int columnIndex = -1;
            if (reader.hasNext()) {
                reader.next();
                for (int c = 0; c < reader.numTokens(); c++) {
                    String value = reader.getString();
                    if ("varCountInSample[align-specific-sample]".equals(value)) {
                        columnIndex = c;
                    }
                }
            }
            while (reader.hasNext()) {
                reader.next();
                for (int c = 0; c < reader.numTokens(); c++) {
                    String value = reader.getString();
                    if (c == columnIndex) {
                        org.junit.Assert.assertEquals("variant count must be 2 in the specific alignment sample", "2", value);

                    }
                }
            }
        }
    }

    @Test
    public void testCleanup() {
        QualityScoreFilter qualityScoreFilter = new QualityScoreFilter();
        LeftOverFilter leftOverFilter = new LeftOverFilter();


        DiscoverVariantPositionData list = new DiscoverVariantPositionData();

        appendInfo(list, 1, (byte) 20, false, 'A', 'C', 0);  // will be filtered by qualityScoreFilter
        appendInfo(list, 1, (byte) 40, false, 'A', 'T', 0);  // will be filtered by leftOverFilter
        appendInfo(list, 10, (byte) 40, true, 'A', 'A', 0);

        ObjectSet<PositionBaseInfo> removed = new ObjectArraySet<PositionBaseInfo>();
        SampleCountInfo[] sampleCounts = sampleCounts(list);
        qualityScoreFilter.filterGenotypes(list, sampleCounts, removed);
        assertEquals(10, sampleCounts[0].counts[SampleCountInfo.BASE_A_INDEX]);
        assertEquals(0, sampleCounts[0].counts[SampleCountInfo.BASE_C_INDEX]);
        assertEquals(1, sampleCounts[0].counts[SampleCountInfo.BASE_T_INDEX]);
        assertEquals(1, removed.size());
        leftOverFilter.filterGenotypes(list, sampleCounts, removed);
        CountFixer fixer = new CountFixer();
        fixer.fix(list, sampleCounts, removed);
        assertEquals(2, removed.size());
        assertEquals(10, list.size());
        assertEquals(10, sampleCounts[0].counts[SampleCountInfo.BASE_A_INDEX]);
        assertEquals(0, sampleCounts[0].counts[SampleCountInfo.BASE_C_INDEX]);
        assertEquals(0, sampleCounts[0].counts[SampleCountInfo.BASE_T_INDEX]);

    }


    @Test
    public void testQuarterFilter() {
        AtLeastAQuarterFilter atLeastAQuarterFilter = new AtLeastAQuarterFilter();


        DiscoverVariantPositionData list = new DiscoverVariantPositionData();

        appendInfo(list, 2, (byte) 40, false, 'A', 'T', 0);  // will be filtered by at least a quarter filter
        appendInfo(list, 3, (byte) 40, false, 'A', 'C', 0);  // will be kept since 3>= 12/4
        appendInfo(list, 4, (byte) 40, false, 'A', 'G', 0);  // will be kept since 4>= 12/4
        appendInfo(list, 12, (byte) 40, true, 'A', 'A', 0);  // allele with most count.

        ObjectSet<PositionBaseInfo> removed = new ObjectArraySet<PositionBaseInfo>();
        SampleCountInfo[] sampleCounts = sampleCounts(list);
        atLeastAQuarterFilter.filterGenotypes(list, sampleCounts, removed);
        assertEquals(2, removed.size());
        CountFixer fixer = new CountFixer();
        fixer.fix(list, sampleCounts, removed);

        assertEquals(19, list.size());
        assertEquals(12, sampleCounts[0].counts[SampleCountInfo.BASE_A_INDEX]);
        assertEquals(3, sampleCounts[0].counts[SampleCountInfo.BASE_C_INDEX]);
        assertEquals(4, sampleCounts[0].counts[SampleCountInfo.BASE_G_INDEX]);
        assertEquals(0, sampleCounts[0].counts[SampleCountInfo.BASE_T_INDEX]);

    }

    @Test
    public void testDiploidFilter() {
        DiploidFilter atLeastAQuarterFilter = new DiploidFilter();


        DiscoverVariantPositionData list = new DiscoverVariantPositionData();

        appendInfo(list, 1, (byte) 40, false, 'A', 'T', 0);  // will be filtered
        appendInfo(list, 4, (byte) 40, false, 'A', 'C', 0);  // will be kept
        appendInfo(list, 2, (byte) 40, false, 'A', 'G', 0);  // will be filtered, two bases
        appendInfo(list, 3, (byte) 40, true, 'A', 'A', 0);  // allele will be kept

        ObjectSet<PositionBaseInfo> removed = new ObjectArraySet<PositionBaseInfo>();
        SampleCountInfo[] sampleCounts = sampleCounts(list);
        atLeastAQuarterFilter.filterGenotypes(list, sampleCounts, removed);
        assertEquals(3, removed.size());
        CountFixer fixer = new CountFixer();
        fixer.fix(list, sampleCounts, removed);

        assertEquals(7, list.size());
        assertEquals(3, sampleCounts[0].counts[SampleCountInfo.BASE_A_INDEX]);
        assertEquals(4, sampleCounts[0].counts[SampleCountInfo.BASE_C_INDEX]);
        assertEquals(0, sampleCounts[0].counts[SampleCountInfo.BASE_G_INDEX]);
        assertEquals(0, sampleCounts[0].counts[SampleCountInfo.BASE_T_INDEX]);

    }

    @Test
    public void testDiploidFilter2() {
        DiploidFilter atLeastAQuarterFilter = new DiploidFilter();


        DiscoverVariantPositionData list = new DiscoverVariantPositionData();

        appendInfo(list, 13, (byte) 40, false, 'A', 'T', 0);  // will be filtered
        appendInfo(list, 12, (byte) 40, false, 'A', 'C', 0);  // will be kept
        appendInfo(list, 13, (byte) 40, false, 'A', 'G', 0);  // will be filtered, two bases
        appendInfo(list, 11, (byte) 40, true, 'A', 'A', 0);  // allele will be kept

        ObjectSet<PositionBaseInfo> removed = new ObjectArraySet<PositionBaseInfo>();
        SampleCountInfo[] sampleCounts = sampleCounts(list);
        atLeastAQuarterFilter.filterGenotypes(list, sampleCounts, removed);
        assertEquals(23, removed.size());
        CountFixer fixer = new CountFixer();
        fixer.fix(list, sampleCounts, removed);

        assertEquals(26, list.size());
        assertEquals(0, sampleCounts[0].counts[SampleCountInfo.BASE_A_INDEX]);
        assertEquals(0, sampleCounts[0].counts[SampleCountInfo.BASE_C_INDEX]);
        assertEquals(13, sampleCounts[0].counts[SampleCountInfo.BASE_G_INDEX]);
        assertEquals(13, sampleCounts[0].counts[SampleCountInfo.BASE_T_INDEX]);

    }

    @Test
    public void testDiploidFilter3() {
        DiploidFilter atLeastAQuarterFilter = new DiploidFilter();


        DiscoverVariantPositionData list = new DiscoverVariantPositionData();

        appendInfo(list, 13, (byte) 40, false, 'A', 'T', 0);  // will be filtered
        appendInfo(list, 12, (byte) 40, false, 'A', 'C', 0);  // will be kept
        appendInfo(list, 11, (byte) 40, false, 'A', 'G', 0);  // will be filtered, two bases
        appendInfo(list, 10, (byte) 40, true, 'A', 'A', 0);  // allele will be kept

        ObjectSet<PositionBaseInfo> removed = new ObjectArraySet<PositionBaseInfo>();
        SampleCountInfo[] sampleCounts = sampleCounts(list);
        atLeastAQuarterFilter.filterGenotypes(list, sampleCounts, removed);
        assertEquals(21, removed.size());
        CountFixer fixer = new CountFixer();
        fixer.fix(list, sampleCounts, removed);

        assertEquals(25, list.size());
        assertEquals(0, sampleCounts[0].counts[SampleCountInfo.BASE_A_INDEX]);
        assertEquals(12, sampleCounts[0].counts[SampleCountInfo.BASE_C_INDEX]);
        assertEquals(0, sampleCounts[0].counts[SampleCountInfo.BASE_G_INDEX]);
        assertEquals(13, sampleCounts[0].counts[SampleCountInfo.BASE_T_INDEX]);

    }

    // construct a SampleCountInfo array holding counts for list
    private SampleCountInfo[] sampleCounts(DiscoverVariantPositionData list) {
        IntSet sampleIndices = new
                IntArraySet();
        for (PositionBaseInfo info : list) {
            sampleIndices.add(info.readerIndex);
        }
        SampleCountInfo[] sci = new SampleCountInfo[sampleIndices.size()];
        for (int i = 0; i < sci.length; i++) {
            sci[i] = new SampleCountInfo();
        }
        for (PositionBaseInfo info : list) {
            final int sampleIndex = info.readerIndex;
            int toBaseIndex = sci[sampleIndex].baseIndex(info.to);
            sci[sampleIndex].counts[toBaseIndex]++;
            if (info.matchesReference) {
                sci[sampleIndex].refCount++;
            } else {
                sci[sampleIndex].varCount++;
            }

        }
        return sci;
    }

    private void appendInfo(DiscoverVariantPositionData list, int n, byte qualityScore, boolean matchesReference, char from, char to, int sampleIndex) {
        for (int i = 0; i < n; i++) {
            PositionBaseInfo info = new PositionBaseInfo();
            info.qualityScore = qualityScore;
            info.from = from;
            info.to = to;
            info.matchesReference = matchesReference;
            info.readerIndex = sampleIndex;
            list.add(info);
        }
    }


    private SampleCountInfo[] makeSampleCounts() {
        SampleCountInfo[] sampleCounts = new SampleCountInfo[1];
        sampleCounts[0] = new SampleCountInfo();

        sampleCounts[0].counts[SampleCountInfo.BASE_A_INDEX] = 5;
        sampleCounts[0].counts[SampleCountInfo.BASE_C_INDEX] = 9;
        sampleCounts[0].counts[SampleCountInfo.BASE_T_INDEX] = 1;
        sampleCounts[0].counts[SampleCountInfo.BASE_OTHER_INDEX] = 1;
        sampleCounts[0].referenceBase = 'A';
        sampleCounts[0].refCount = 5;
        sampleCounts[0].varCount = 11;
        return sampleCounts;
    }

    private SampleCountInfo[] makeTwoSampleCounts() {
        SampleCountInfo[] sampleCounts = new SampleCountInfo[2];
        sampleCounts[0] = new SampleCountInfo();
        sampleCounts[0].counts[SampleCountInfo.BASE_A_INDEX] = 5;
        sampleCounts[0].counts[SampleCountInfo.BASE_T_INDEX] = 1;
        sampleCounts[0].counts[SampleCountInfo.BASE_C_INDEX] = 9;
        sampleCounts[0].counts[SampleCountInfo.BASE_OTHER_INDEX] = 1;
        sampleCounts[0].referenceBase = 'A';
        sampleCounts[0].refCount = 5;
        sampleCounts[0].varCount = 11;
        sampleCounts[0].sampleIndex = 0;

        sampleCounts[1] = new SampleCountInfo();
        sampleCounts[1].counts[SampleCountInfo.BASE_A_INDEX] = 10;
        sampleCounts[1].counts[SampleCountInfo.BASE_T_INDEX] = 4;
        sampleCounts[1].counts[SampleCountInfo.BASE_C_INDEX] = 0;
        sampleCounts[1].counts[SampleCountInfo.BASE_OTHER_INDEX] = 2;
        sampleCounts[1].referenceBase = 'A';
        sampleCounts[1].refCount = 10;
        sampleCounts[1].varCount = 6;
        sampleCounts[1].sampleIndex = 1;
        return sampleCounts;
    }


    @Test
    public void testGenotypes() throws IOException, JSAPException {
        DiscoverSequenceVariantsMode mode = new DiscoverSequenceVariantsMode();
        int i = 1;
        String outputFilename = "out-genotypes-" + i + ".tsv";
        String[] args = constructArgumentString(
                basenames, BASE_TEST_DIR + "/" + outputFilename, "samples").split("[\\s]");
        args = add(args, new String[]{"--format", DiscoverSequenceVariantsMode.OutputFormat.GENOTYPES.toString()});

        mode.setDisableAtLeastQuarterFilter(true);
        configureTestGenome(mode);
        mode.configure(args);

        mode.execute();
        assertEquals(new File("test-data/discover-variants/expected-output-genotypes.tsv"),
                new File(BASE_TEST_DIR + "/" + outputFilename)
        );

    }

    @Test
    public void testQualityScoreAdjuster() {
        QualityScoreFilter adjuster = new QualityScoreFilter();
        SampleCountInfo[] sampleCounts = makeSampleCounts();
        IntArrayList scores = IntArrayList.wrap(new int[]{10, 20, 30, 40, 40, 40, 40});

        final DiscoverVariantPositionData list = makeListWithScores(sampleCounts, scores);
        assertEquals(16, list.size());
        ObjectSet<PositionBaseInfo> filteredList = new ObjectArraySet<PositionBaseInfo>();
        adjuster.filterGenotypes(list, sampleCounts, filteredList);

        System.out.println("list: " + list);
        assertEquals(16, list.size());
        assertEquals(4, filteredList.size());
        assertEquals(5, sampleCounts[0].refCount);
        assertEquals(7, sampleCounts[0].varCount);

    }

    @Test
    public void testAdjustVarCount() {
        QualityScoreFilter adjuster1 = new QualityScoreFilter();
        LeftOverFilter adjuster2 = new LeftOverFilter();
        SampleCountInfo[] sampleCounts = makeTwoSampleCounts();

        assertEquals(5, sampleCounts[0].refCount);
        assertEquals(11, sampleCounts[0].varCount);
        assertEquals(10, sampleCounts[1].refCount);
        assertEquals(6, sampleCounts[1].varCount);

        IntArrayList scores = IntArrayList.wrap(new int[]{10, 20, 30, 40, 40, 40, 40});

        final DiscoverVariantPositionData list = makeListWithScores(sampleCounts, scores);
        assertEquals(32, list.size());
        ObjectSet<PositionBaseInfo> filteredList = new ObjectArraySet<PositionBaseInfo>();
        adjuster1.filterGenotypes(list, sampleCounts, filteredList);
        adjuster2.filterGenotypes(list, sampleCounts, filteredList);

        System.out.println("list: " + list);
        System.out.println("filtered: " + filteredList);

        assertEquals(22, filteredList.size());
        assertEquals(0, sampleCounts[0].refCount);
        assertEquals(0, sampleCounts[0].varCount);

        assertEquals(10, sampleCounts[1].refCount);
        assertEquals(6 - 6, sampleCounts[1].varCount);

    }

    @Test
    public void testAdjustVarCount2() {
        QualityScoreFilter adjuster1 = new QualityScoreFilter();
        LeftOverFilter adjuster2 = new LeftOverFilter();
        SampleCountInfo[] sampleCounts = makeTwoSampleCounts();

        assertEquals(5, sampleCounts[0].refCount);
        assertEquals(11, sampleCounts[0].varCount);

        assertEquals(10, sampleCounts[1].refCount);
        assertEquals(6, sampleCounts[1].varCount);

        IntArrayList scores = IntArrayList.wrap(new int[]{40, 40, 40, 40, 40, 40, 10, 40, 40, 40, 40, 40, 40, 40, 40, 40});

        final DiscoverVariantPositionData list = makeListWithScores(sampleCounts, scores);
        assertEquals(32, list.size());
        ObjectSet<PositionBaseInfo> filteredList = new ObjectArraySet<PositionBaseInfo>();
        adjuster1.filterGenotypes(list, sampleCounts, filteredList);
        adjuster2.filterGenotypes(list, sampleCounts, filteredList);

        System.out.println("list: " + list);
        System.out.println("filtered: " + filteredList);

        assertEquals(3, filteredList.size());
        assertEquals(5, sampleCounts[0].refCount);
        assertEquals(8, sampleCounts[0].varCount);

        assertEquals(10, sampleCounts[1].refCount);
        assertEquals(6, sampleCounts[1].varCount);

    }

    private ObjectArrayList<ReadIndexStats> makeReadIndexStats() {
        ObjectArrayList<ReadIndexStats> result = new ObjectArrayList<ReadIndexStats>();
        ReadIndexStats stat = new ReadIndexStats();

        stat.basename = "basename1";
        stat.readerIndex = 0;
        stat.countReferenceBases = new long[]{1000, 1000, 1000};
        stat.countVariationBases = new long[]{100, 500, 900};
        result.add(stat);

        stat = new ReadIndexStats();
        stat.basename = "basename2";
        stat.readerIndex = 0;
        stat.countReferenceBases = new long[]{1000, 1000, 1000};
        stat.countVariationBases = new long[]{10, 10, 10};
        result.add(stat);
        return result;
    }


    private DiscoverVariantPositionData makeList(SampleCountInfo[] sampleCounts, IntArrayList readIndices) {
        IntIterator nextReadIndexIterator = readIndices.iterator();
        DiscoverVariantPositionData list = new DiscoverVariantPositionData();
        for (SampleCountInfo sampleInfo : sampleCounts) {
            for (int baseIndex = 0; baseIndex < SampleCountInfo.BASE_MAX_INDEX; baseIndex++) {

                for (int i = 0; i < sampleInfo.counts[baseIndex]; i++) {

                    PositionBaseInfo info = new PositionBaseInfo();
                    final char base = sampleInfo.base(baseIndex);
                    info.to = base;
                    if (base == 'A') {
                        info.matchesReference = true;
                        info.from = base;
                        info.matchesForwardStrand = true;
                    }
                    info.readerIndex = sampleInfo.sampleIndex;
                    if (!nextReadIndexIterator.hasNext()) {

                        // wrap back to the start of read indices:
                        nextReadIndexIterator = readIndices.iterator();
                    }
                    info.readIndex = nextReadIndexIterator.nextInt();
                    list.add(info);
                    System.out.println("info: " + info);
                }

            }
        }
        return list;
    }

    private DiscoverVariantPositionData makeListWithScores(SampleCountInfo[] sampleCounts, IntArrayList qualityScores) {
        IntIterator nextQualityIterator = qualityScores.iterator();
        DiscoverVariantPositionData list = new DiscoverVariantPositionData();
        for (SampleCountInfo sampleInfo : sampleCounts) {
            for (int baseIndex = 0; baseIndex < SampleCountInfo.BASE_MAX_INDEX; baseIndex++) {

                for (int i = 0; i < sampleInfo.counts[baseIndex]; i++) {

                    PositionBaseInfo info = new PositionBaseInfo();

                    final char base = sampleInfo.base(baseIndex);
                    info.to = base;
                    if (base == 'A') {
                        info.matchesReference = true;
                        info.from = base;
                        info.matchesForwardStrand = true;
                    }
                    info.readerIndex = sampleInfo.sampleIndex;
                    if (!nextQualityIterator.hasNext()) {

                        // wrap back to the start of read indices:
                        nextQualityIterator = qualityScores.iterator();
                    }
                    info.readIndex = 1;
                    info.qualityScore = (byte) nextQualityIterator.nextInt();
                    list.add(info);
                    System.out.println("info: " + info);
                }

            }
        }
        return list;
    }


    private int nextReadIndex() {
        return 0;  //To change body of created methods use File | Settings | File Templates.
    }

    private String[] shuffle
            (String[] strings) {

        final ObjectArrayList<String> list = ObjectArrayList.wrap(strings);
        Collections.shuffle(list);
        return list.toArray(new String[list.size()]);
    }

    private String[] add
            (String[] basenames, String[] specificAlignments) {
        String result[] = new String[basenames.length + specificAlignments.length];
        int i = 0;
        for (String s : basenames) {
            result[i++] = s;
        }
        for (String s : specificAlignments) {
            result[i++] = s;
        }
        return result;
    }

    @AfterClass
    public static void cleanupTestDirectory
            () throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Deleting base test directory: " + BASE_TEST_DIR);
        }
        //    FileUtils.forceDeleteOnExit(new File(BASE_TEST_DIR));
    }

    @Before
    public void setUp
            () throws IOException {

        final File dir = new File(BASE_TEST_DIR);
        dir.mkdirs();


        final ReadsWriter referenceWriter = new ReadsWriter(FileUtils.openOutputStream(
                new File("test-results/alignments/last-to-compact/last-reference.compact-reads")));
        referenceWriter.setIdentifier("0");
        referenceWriter.appendEntry();
        referenceWriter.close();


    }

    private static String constructArgumentString(String[] basenames, Object outputfilename, Object evalString) {
        MutableString basenamesString = new MutableString();
        for (String basename : basenames) {
            basenamesString.append(basename).append(" ");
        }

        MutableString groups = new MutableString();
        groups.append("A=");
        for (int i = 0; i < NUM_BASENAME_PER_GROUP; i++) {
            groups.append(FilenameUtils.getBaseName(basenames[i]));
            if (i != NUM_BASENAME_PER_GROUP - 1) {
                groups.append(",");
            }
        }
        groups.append("/B=");
        for (int i = NUM_BASENAME_PER_GROUP; i < NUM_BASENAME_PER_GROUP * 2; i++) {
            groups.append(FilenameUtils.getBaseName(basenames[i]));
            if (i != NUM_BASENAME_PER_GROUP * 2 - 1) {
                groups.append(",");
            }
        }

        return String.format("--mode discover-sequence-variants " +
                "--variation-stats %s " +
                "--groups %s " +
                "--compare A/B " +
                "--genome use-dummy-in-test " +
                "--eval %s " +
                "--minimum-variation-support 1 " +
                "--threshold-distinct-read-indices 1 " +
                "--output %s " +
                "%s", statsFilename, groups, evalString, outputfilename, basenamesString);
    }

    @BeforeClass
    public static void initializeTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Creating base test directory: " + BASE_TEST_DIR);
        }

        FileUtils.forceMkdir(new File(BASE_TEST_DIR));
        PrintWriter outTSV = new PrintWriter(statsFilename);
        outTSV.println("basename\tread-index\tcount-variation-bases\tbases-at-index/all-variations-bases\tbases-at-index/all-reference-bases\tcount-reference-bases\tcount-reference-bases-at-index");
        for (int i = 0; i < NUM_BASENAME_PER_GROUP * 2; i++) {
            basenames[i] = BASE_TEST_DIR + "/basen" + i;
            writeAlignment(basenames, i, 'G');
            appendToStats(basenames, i, outTSV);
        }
        specificAlignments = new String[1];
        // write a new alignment that is different from all others:
        specificAlignments[0] = BASE_TEST_DIR + "/align-specific-sample";
        writeAlignment(specificAlignments, 0, 'C');
        appendToStats(specificAlignments, 0, outTSV);

        outTSV.close();

    }

    private static void appendToStats(String[] basenames, int basenameIndex, PrintWriter outTSV) throws IOException {
        String tokens = "basename\tread-index\tcount-variation-bases\tbases-at-index/all-variations-bases\tbases-at-index/all-reference-bases\tcount-reference-bases\tcount-reference-bases-at-index\n" +
                "basename0\t1\t2924830\t0.097711641\t-0.004026\t-726568406\t68804818\n" +
                "basename0\t2\t1683737\t0.056249664\t-0.002317\t-726568406\t70050107\n" +
                "basename0\t3\t973088\t0.032508565\t-0.001339\t-726568406\t70793744\n" +
                "basename0\t4\t624546\t0.020864602\t-0.00086\t-726568406\t71124447\n" +
                "basename0\t5\t458581\t0.015320105\t-0.000631\t-726568406\t71275227\n" +
                "basename0\t6\t455194\t0.015206953\t-0.000626\t-726568406\t71268937\n" +
                "basename0\t7\t338682\t0.011314563\t-0.000466\t-726568406\t71383137\n" +
                "basename0\t8\t342731\t0.011449831\t-0.000472\t-726568406\t71377688\n" +
                "basename0\t9\t305846\t0.01021759\t-0.000421\t-726568406\t71414272\n" +
                "basename0\t10\t391497\t0.013078987\t-0.000539\t-726568406\t71328565\n" +
                "basename0\t11\t285087\t0.009524081\t-0.000392\t-726568406\t71434921\n" +
                "basename0\t12\t448654\t0.014988467\t-0.000617\t-726568406\t71271330\n" +
                "basename0\t13\t263024\t0.008787009\t-0.000362\t-726568406\t71456737\n" +
                "basename0\t14\t259276\t0.008661797\t-0.000357\t-726568406\t71460537\n" +
                "basename0\t15\t255673\t0.008541429\t-0.000352\t-726568406\t71463985\n" +
                "basename0\t16\t265280\t0.008862376\t-0.000365\t-726568406\t71454330\n" +
                "basename0\t17\t286275\t0.00956377\t-0.000394\t-726568406\t71433717\n" +
                "basename0\t18\t253905\t0.008482364\t-0.000349\t-726568406\t71466066\n" +
                "basename0\t19\t242706\t0.008108232\t-0.000334\t-726568406\t71477255\n" +
                "basename0\t20\t255099\t0.008522253\t-0.000351\t-726568406\t71464680\n" +
                "basename0\t21\t244071\t0.008153834\t-0.000336\t-726568406\t71475793\n" +
                "basename0\t22\t243361\t0.008130114\t-0.000335\t-726568406\t71476488\n" +
                "basename0\t23\t237500\t0.007934312\t-0.000327\t-726568406\t71482128\n" +
                "basename0\t24\t393818\t0.013156526\t-0.000542\t-726568406\t71325748\n" +
                "basename0\t25\t290470\t0.009703915\t-0.0004\t-726568406\t71429370\n" +
                "basename0\t26\t283203\t0.009461141\t-0.00039\t-726568406\t71436438\n" +
                "basename0\t27\t398203\t0.013303019\t-0.000548\t-726568406\t71321346\n" +
                "basename0\t28\t244253\t0.008159914\t-0.000336\t-726568406\t71475340\n" +
                "basename0\t29\t250658\t0.00837389\t-0.000345\t-726568406\t71469218\n" +
                "basename0\t30\t252561\t0.008437465\t-0.000348\t-726568406\t71466726\n" +
                "basename0\t31\t270395\t0.009033256\t-0.000372\t-726568406\t71449213\n" +
                "basename0\t32\t255815\t0.008546173\t-0.000352\t-726568406\t71463771\n" +
                "basename0\t33\t264515\t0.008836819\t-0.000364\t-726568406\t71455248\n" +
                "basename0\t34\t300673\t0.010044773\t-0.000414\t-726568406\t71419253\n" +
                "basename0\t35\t297241\t0.009930118\t-0.000409\t-726568406\t71422416\n" +
                "basename0\t36\t295019\t0.009855886\t-0.000406\t-726568406\t71424664\n" +
                "basename0\t37\t307167\t0.010261722\t-0.000423\t-726568406\t71412201\n" +
                "basename0\t38\t310107\t0.01035994\t-0.000427\t-726568406\t71409148\n" +
                "basename0\t39\t491224\t0.01641063\t-0.000676\t-726568406\t71228249\n" +
                "basename0\t40\t339371\t0.011337581\t-0.000467\t-726568406\t71380591\n" +
                "basename0\t41\t447670\t0.014955594\t-0.000616\t-726568406\t71273454\n" +
                "basename0\t42\t381186\t0.012734521\t-0.000525\t-726568406\t71343713\n" +
                "basename0\t43\t445198\t0.01487301\t-0.000613\t-726568406\t71287803\n" +
                "basename0\t44\t482766\t0.016128068\t-0.000664\t-726568406\t71264200\n" +
                "basename0\t45\t656445\t0.021930272\t-0.000903\t-726568406\t71075738\n" +
                "basename0\t46\t704999\t0.023552346\t-0.00097\t-726568406\t71016151\n" +
                "basename0\t47\t1049645\t0.035066153\t-0.001445\t-726568406\t70666569\n" +
                "basename0\t48\t1674529\t0.055942047\t-0.002305\t-726568406\t70040028\n" +
                "basename0\t49\t2453779\t0.081974943\t-0.003377\t-726568406\t69260573\n" +
                "basename0\t50\t3343931\t0.111712812\t-0.004602\t-726568406\t68370412\n" +
                "basename0\t51\t1009797\t0.033734925\t-0.00139\t-726568406\t11172400\n";
        StringReader reader = new StringReader(tokens);
        BufferedReader br = new BufferedReader(reader);
        String line;

        // skip header line:
        br.readLine();

        while ((line = br.readLine()) != null) {
            String[] toks = line.split("[\\t]");

            outTSV.print(FilenameUtils.getBaseName(basenames[basenameIndex]));

            for (int i = 1; i < toks.length; i++) {

                outTSV.print("\t");
                outTSV.print(toks[i]);

            }
            outTSV.println();

        }

    }

    private static void writeAlignment(String[] basenames, int basenameIndex, char toBase) throws IOException {

        AlignmentWriter writer = new AlignmentWriterImpl(basenames[basenameIndex]);
        writer.setTargetLengths(new int[]{10000});
        writer.setSorted(true);
        writer.setTargetIdentifiersArray(new String[]{"target1"});

        int numAlignmentEntries = 10;
        int readIndex = 0;
        int referencePosition = 1;
        int positionStart = 100;
        writeAlignmentEntries(toBase, writer, numAlignmentEntries, readIndex, referencePosition, positionStart);
        if (toBase == 'C') {
            referencePosition = 1;
            writeAlignmentEntries(toBase, writer, numAlignmentEntries, readIndex, referencePosition, positionStart);
        }
        writer.close();
    }

    private static void writeAlignmentEntries(char toBase, AlignmentWriter writer, int numAlignmentEntries, int readIndex, int referencePosition, int positionStart) throws IOException {
        for (int j = 0; j < numAlignmentEntries; j++) {
            final Alignments.AlignmentEntry.Builder builder =
                    Alignments.AlignmentEntry.newBuilder();

            builder.setQueryIndex(readIndex++);
            builder.setTargetIndex(0);

            builder.setQueryLength(50);
            builder.setPosition(referencePosition++ + positionStart);
            builder.setMatchingReverseStrand(j % 2 == 0);
            //  builder.setMatchingReverseStrand(false);

            builder.setScore(50);
            builder.setNumberOfIndels(0);

            builder.setQueryAlignedLength(50);
            int multiplicity = 1;
            builder.setMultiplicity(multiplicity);
            builder.setTargetAlignedLength(50);
            Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder();
            varBuilder.setFrom("A");
            varBuilder.setTo(Character.toString(toBase));

            varBuilder.setToQuality(ByteString.copyFrom(new byte[]{40}));
            varBuilder.setPosition(25);
            varBuilder.setReadIndex(25);

            //   System.out.printf("%s j=%d var A/G at %d%n", basenames[basenameIndex], j, 25 + referencePosition + positionStart - 1);
            builder.addSequenceVariations(varBuilder);
            final Alignments.AlignmentEntry entry = builder.build();
            writer.appendEntry(entry);

        }
    }


}
