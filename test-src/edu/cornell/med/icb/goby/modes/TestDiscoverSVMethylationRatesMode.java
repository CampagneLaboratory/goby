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
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import edu.cornell.med.icb.goby.reads.ReadsWriterImpl;
import edu.cornell.med.icb.goby.util.TestFiles;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
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
import java.util.ArrayList;
import java.util.Collections;

import static org.junit.Assert.assertTrue;

/**
 * @author Fabien Campagne
 *         Date: Mar 12, 2011
 *         Time: 12:15:11 PM
 */
public class TestDiscoverSVMethylationRatesMode extends TestFiles {


    private static final Log LOG = LogFactory.getLog(TestDiscoverSVMethylationRatesMode.class);
    private static final String BASE_TEST_DIR = "test-results/discover-variants-methylation";
    private static final int NUM_BASENAME_PER_GROUP = 5;
    private static final String[] basenames = new String[NUM_BASENAME_PER_GROUP * 2];
    private static final int NUM_TRIALS = 10;
    private static String statsFilename = BASE_TEST_DIR + "/variations-stats2.tsv";
    private static String[] specificAlignments;
    private String header = "##fileformat=VCFv4.1\n" +
            "##Goby=UNKNOWN\n" +
            "##INFO=<ID=BIOMART_COORDS,Number=1,Type=String,Description=\"Coordinates for use with Biomart.\">\n" +
            "##INFO=<ID=Strand,Number=1,Type=String,Description=\"Strand of the cytosine site on the reference sequence.\">\n" +
            "##INFO=<ID=LOD[methylated/not-so],Number=1,Type=Float,Description=\"Log2 of the odds-ratio of observing methylation in  group methylated versus group not-so\">\n" +
            "##INFO=<ID=LOD_SE[methylated/not-so],Number=1,Type=Float,Description=\"Standard Error of the log2 of the odds-ratio between group methylated and group not-so\">\n" +
            "##INFO=<ID=LOD_Z[methylated/not-so],Number=1,Type=Float,Description=\"Z value of the odds-ratio between group methylated and group not-so\">\n" +
            "##INFO=<ID=FisherP[methylated/not-so],Number=1,Type=Float,Description=\"Fisher exact P-value of observing as large a difference by chance between group methylated and group not-so.\">\n" +
            "##INFO=<ID=#C Group[methylated],Number=1,Type=Integer,Description=\"Number of unmethylated Cytosines at site in group methylated.\">\n" +
            "##INFO=<ID=#Cm Group[methylated],Number=1,Type=Integer,Description=\"Number of methylated Cytosines at site in group methylated.\">\n" +
            "##INFO=<ID=#C Group[not-so],Number=1,Type=Integer,Description=\"Number of unmethylated Cytosines at site in group not-so.\">\n" +
            "##INFO=<ID=#Cm Group[not-so],Number=1,Type=Integer,Description=\"Number of methylated Cytosines at site in group not-so.\">\n" +
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth of sequencing across groups at this site\">\n" +
            "##FORMAT=<ID=MR,Number=1,Type=Integer,Description=\"Methylation rate. 0-100%, 100% indicate fully methylated.\">\n" +
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
            "##FORMAT=<ID=BC,Number=1,Type=String,Description=\"Base counts in format A=?;T=?;C=?;G=?;N=?.\">\n" +
            "##FORMAT=<ID=GB,Number=1,Type=String,Description=\"Number of bases that pass base filters in this sample.\">\n" +
            "##FORMAT=<ID=FB,Number=1,Type=String,Description=\"Number of bases that failed base filters in this sample.\">\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmethylated\tnot-so\n";

    RandomAccessSequenceTestSupport genome;

    private String[] sequences = {
            "ATTTACCGG", //CpG
            "TAGATACGGAT",//CpG
            "ACTCTAGACTA", //CpT
            "CATTTTGCAAC", //CpA
            "ATCTATGCCTA", //CpG - negative
    };



    @Test
    public void testMethylationFormatSwap() throws IOException, JSAPException {
        final MethylationRateVCFOutputFormat outputFormat = new MethylationRateVCFOutputFormat();
        DiscoverVariantIterateSortedAlignments iterator = new DiscoverVariantIterateSortedAlignments(outputFormat) {
            @Override
            public CharSequence getReferenceId(int targetIndex) {
                return "ref-id";
            }
        };
        DiscoverSequenceVariantsMode mode = new DiscoverSequenceVariantsMode() {
            @Override
            public String[] getSamples() {
                return new String[]{"not-so", "methylated"};
            }

            @Override
            public String[] getGroups() {
                return new String[]{"not-so", "methylated",};
            }

            @Override
            public int[] getReaderIndexToGroupIndex() {
                return new int[]{0, 1}; // sample index is group index;
            }

            @Override
            public ArrayList<GroupComparison> getGroupComparisons() {
                ArrayList<GroupComparison> list = new ArrayList<GroupComparison>();
                list.add(new GroupComparison("methylated", "not-so", 0, 1, 0));
                return list;
            }
        };
        StringWriter writer = new StringWriter();
        outputFormat.allocateStorage(2, 2);
        outputFormat.defineColumns(new PrintWriter(writer), mode);
        outputFormat.setGenome(genome);

        outputFormat.setMinimumEventThreshold(0);
        outputFormat.writeRecord(iterator, makeTwoSampleCounts(), 0, 0, list4(), 0, 1);
        String stringB = writer.getBuffer().toString();
        assertTrue(stringB, stringB.contains("1/2:A=10,T=4,N=2:16:0:33"));
        assertTrue(stringB, stringB.contains("0/1/2:A=5,T=1,C=9,N=1:16:0:0"));
        assertTrue(stringB, stringB.contains("#Cm_Group[methylated]=10;"));
        assertTrue(stringB, stringB.contains("#C_Group[methylated]=20;"));
        // check biomart-span included in result:
        assertTrue("biomart span must be included in result: " + stringB, stringB.contains("ref-id:1:1"));
    }

    @Test
    public void testMethylationFormat() throws IOException, JSAPException {
        final MethylationRateVCFOutputFormat outputFormat = new MethylationRateVCFOutputFormat();
        DiscoverVariantIterateSortedAlignments iterator = new DiscoverVariantIterateSortedAlignments(outputFormat) {
            @Override
            public CharSequence getReferenceId(int targetIndex) {
                return "ref-id";
            }
        };
        DiscoverSequenceVariantsMode mode = new DiscoverSequenceVariantsMode() {
            @Override
            public String[] getSamples() {
                return new String[]{"methylated", "not-so"};
            }

            @Override
            public String[] getGroups() {
                return new String[]{"methylated", "not-so"};
            }

            @Override
            public int[] getReaderIndexToGroupIndex() {
                return new int[]{0, 1}; // sample index is group index;
            }

            @Override
            public ArrayList<GroupComparison> getGroupComparisons() {
                ArrayList<GroupComparison> list = new ArrayList<GroupComparison>();
                list.add(new GroupComparison("methylated", "not-so", 0, 1, 0));
                return list;
            }
        };
        StringWriter writer = new StringWriter();
        outputFormat.allocateStorage(2, 2);
        outputFormat.setGenome(genome);
        outputFormat.defineColumns(new PrintWriter(writer), mode);
        outputFormat.setMinimumEventThreshold(0);
        outputFormat.writeRecord(iterator, makeTwoSampleCounts(), 0, 0, list1(), 0, 1);
        writer.flush();

        String stringB = writer.getBuffer().toString();
        assertTrue(stringB, stringB.contains("#Cm_Group[methylated]=1;"));
        assertTrue(stringB, stringB.contains("#Cm_Group[not-so]=0;"));


        writer = new StringWriter();
        outputFormat.allocateStorage(2, 2);
        outputFormat.defineColumns(new PrintWriter(writer), mode);

        final SampleCountInfo[] sampleCounts = makeTwoSampleCounts();
        sampleCounts[0].referenceBase = 'G';
        outputFormat.writeRecord(iterator, sampleCounts, 0, 0, list2(), 0, 1);
        writer.flush();
        stringB = writer.getBuffer().toString();
        assertTrue(stringB, stringB.contains("#Cm_Group[methylated]=1;"));
        assertTrue(stringB, stringB.contains("#Cm_Group[not-so]=0;"));
        assertTrue(stringB, stringB.contains("MR[methylated]=100.0;"));
        assertTrue(stringB, stringB.contains("MR[not-so]=NaN;"));


        writer = new StringWriter();

        outputFormat.allocateStorage(2, 2);
        outputFormat.defineColumns(new PrintWriter(writer), mode);

        outputFormat.writeRecord(iterator, makeTwoSampleCounts(), 0, 0, list3(), 0, 1);
        stringB = writer.getBuffer().toString();
        assertTrue(stringB, stringB.contains("#Cm_Group[methylated]=1;"));
        assertTrue(stringB, stringB.contains("#C_Group[methylated]=1;"));
        assertTrue(stringB, stringB.contains("#Cm_Group[not-so]=0;"));
        assertTrue(stringB, stringB.contains("#C_Group[not-so]=0;"));


        /*  org.junit.Assert.assertEquals("content must match", header+
                "ref-id\t1\t\tC\tA,T,N\t\t\tBIOMART_COORDS=ref-id:1:1;Strand=+;LOD[methylated/not-so]=NaN;LOD_SE[methylated/not-so]=NaN;LOD_Z[methylated/not-so]=NaN;FisherP[methylated/not-so]=NaN;#C Group[methylated]=0;#Cm Group[methylated]=1;#C Group[not-so]=0;#Cm Group[not-so]=0;DP=1\tMR:GT:BC:GB:FB\t100:0/1/2/3:A=5,T=1,C=9,G=0,N=1:16:0\t0:1/2/3:A=10,T=4,C=0,G=0,N=2:16:0",
                buffer.toString());
        */
    }

    @Test
    public void testMethylationFormatGenomicContext() throws IOException, JSAPException {
        final MethylationRateVCFOutputFormat outputFormat = new MethylationRateVCFOutputFormat();

        DiscoverVariantIterateSortedAlignments iterator = new DiscoverVariantIterateSortedAlignments(outputFormat) {
            @Override
            public CharSequence getReferenceId(int targetIndex) {
               return Integer.toString(targetIndex);
            }
        };
        DiscoverSequenceVariantsMode mode = new DiscoverSequenceVariantsMode() {
            @Override
            public String[] getSamples() {
                return new String[]{"methylated", "not-so"};
            }

            @Override
            public String[] getGroups() {
                return new String[]{"methylated", "not-so"};
            }

            @Override
            public int[] getReaderIndexToGroupIndex() {
                return new int[]{0, 1}; // sample index is group index;
            }

            @Override
            public ArrayList<GroupComparison> getGroupComparisons() {
                ArrayList<GroupComparison> list = new ArrayList<GroupComparison>();
                list.add(new GroupComparison("methylated", "not-so", 0, 1, 0));
                return list;
            }
        };
        outputFormat.setGenome(genome);

        StringWriter writer = new StringWriter();
        outputFormat.allocateStorage(2, 2);
        outputFormat.defineColumns(new PrintWriter(writer), mode);
        outputFormat.setMinimumEventThreshold(0);
        outputFormat.writeRecord(iterator, makeTwoSampleCounts(), 1,6, list1(), 0, 1);
        writer.flush();

        String stringB = writer.getBuffer().toString();
        assertTrue(stringB, stringB.contains("Context=CpG"));

        writer = new StringWriter();
        outputFormat.allocateStorage(2, 2);
        outputFormat.defineColumns(new PrintWriter(writer), mode);
        outputFormat.setMinimumEventThreshold(0);

        final int referenceIndex = 2;
        outputFormat.setGenomeReferenceIndex(referenceIndex);
        outputFormat.writeRecord(iterator, makeTwoSampleCounts(), referenceIndex, 1, list1(), 0, 1);
        writer.flush();

        String stringC = writer.getBuffer().toString();
        assertTrue(stringC , stringC.contains("Context=CpT"));


        writer = new StringWriter();
        outputFormat.allocateStorage(2, 2);
        outputFormat.defineColumns(new PrintWriter(writer), mode);
        outputFormat.setMinimumEventThreshold(0);
        final int referenceIndex1 = 1;
        outputFormat.setGenomeReferenceIndex(referenceIndex1);
        outputFormat.writeRecord(iterator, makeTwoSampleCounts(), referenceIndex1, 2, list1(), 0, 1);
        stringB = writer.getBuffer().toString();
        assertTrue(stringB, stringB.contains("CpT"));

        writer = new StringWriter();
        outputFormat.allocateStorage(2, 2);
        outputFormat.defineColumns(new PrintWriter(writer), mode);
        outputFormat.setMinimumEventThreshold(0);
        final int referenceIndex4 = 4;
        outputFormat.setGenomeReferenceIndex(referenceIndex4);
        outputFormat.writeRecord(iterator, makeTwoSampleCounts(), referenceIndex4, 7, list1(), 0, 1);
        stringB = writer.getBuffer().toString();
        assertTrue(stringB, stringB.contains("CpC"));


    }

    private DiscoverVariantPositionData list4() {
        DiscoverVariantPositionData result = new DiscoverVariantPositionData();
        PositionBaseInfo info;
        // make 33% methylation in sampleIndex=1 "methylated"
        for (int i = 0; i < 10; i++) {
            info = new PositionBaseInfo();
            info.matchesReference = true;
            info.matchesForwardStrand = true;
            info.from = 'C';
            info.to = 'C';
            info.readerIndex = 1;
            result.add(info);
        }
        for (int i = 0; i < 20; i++) {
            info = new PositionBaseInfo();
            info.matchesReference = false;      // not methylated
            info.matchesForwardStrand = true;
            info.from = 'C';
            info.to = 'T';               // not methylated
            info.readerIndex = 1;
            result.add(info);
        }
        return result;
    }

    private DiscoverVariantPositionData list3() {
        DiscoverVariantPositionData result = new DiscoverVariantPositionData();
        PositionBaseInfo info = new PositionBaseInfo();
        info.matchesReference = true;
        info.matchesForwardStrand = true;
        info.from = 'C';
        info.to = 'C';

        result.add(info);

        info = new PositionBaseInfo();
        info.matchesReference = false;
        info.matchesForwardStrand = true;
        info.from = 'C';
        info.to = 'T';

        result.add(info);

        return result;
    }

    private DiscoverVariantPositionData list2() {
        DiscoverVariantPositionData result = new DiscoverVariantPositionData();
        PositionBaseInfo info = new PositionBaseInfo();
        info.matchesReference = true;
        info.matchesForwardStrand = false;
        info.from = 'G';
        info.to = 'G';

        result.add(info);


        return result;
    }

    private DiscoverVariantPositionData list1() {
        DiscoverVariantPositionData result = new DiscoverVariantPositionData();
        final PositionBaseInfo info = new PositionBaseInfo();
        info.matchesReference = true;
        info.matchesForwardStrand = true;
        info.from = 'C';
        info.to = 'C';

        result.add(info);
        return result;
    }


    private SampleCountInfo[] makeTwoSampleCounts() {
        SampleCountInfo[] sampleCounts = new SampleCountInfo[2];
        sampleCounts[0] = new SampleCountInfo();
        sampleCounts[0].setGenotypeCount(SampleCountInfo.BASE_A_INDEX, 5);
        sampleCounts[0].setGenotypeCount(SampleCountInfo.BASE_T_INDEX,1);
        sampleCounts[0].setGenotypeCount(SampleCountInfo.BASE_C_INDEX,9);
        sampleCounts[0].setGenotypeCount(SampleCountInfo.BASE_OTHER_INDEX,1);
        sampleCounts[0].referenceBase = 'C';
        sampleCounts[0].refCount = 5;
        sampleCounts[0].varCount = 11;
        sampleCounts[0].sampleIndex = 0;

        sampleCounts[1] = new SampleCountInfo();
        sampleCounts[1].setGenotypeCount(SampleCountInfo.BASE_A_INDEX,10);
        sampleCounts[1].setGenotypeCount(SampleCountInfo.BASE_T_INDEX,4);
        sampleCounts[1].setGenotypeCount(SampleCountInfo.BASE_C_INDEX,0);
        sampleCounts[1].setGenotypeCount(SampleCountInfo.BASE_OTHER_INDEX, 2);
        sampleCounts[1].referenceBase = 'C';
        sampleCounts[1].refCount = 10;
        sampleCounts[1].varCount = 6;
        sampleCounts[1].sampleIndex = 1;
        return sampleCounts;
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

                for (int i = 0; i < sampleInfo.getGenotypeCount(baseIndex); i++) {

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

         genome = new RandomAccessSequenceTestSupport(sequences);
        final ReadsWriter referenceWriter = new ReadsWriterImpl(FileUtils.openOutputStream(
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
                "--eval %s " +
                "--vcf " +
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