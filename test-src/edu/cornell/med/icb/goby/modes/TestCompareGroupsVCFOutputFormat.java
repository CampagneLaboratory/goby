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

import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.goby.Release1_9_7_2;
import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import edu.cornell.med.icb.goby.alignments.DiscoverVariantIterateSortedAlignments;
import edu.cornell.med.icb.goby.alignments.PositionBaseInfo;
import edu.cornell.med.icb.goby.alignments.SampleCountInfo;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import org.easymock.EasyMock;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.PrintWriter;
import java.io.StringWriter;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;
import static org.easymock.EasyMock.*;

/**
 * @author Fabien Campagne
 *         Date: 5/30/11
 *         Time: 10:37 AM
 */
@SuppressWarnings({"IOResourceOpenedButNotSafelyClosed"})
public class TestCompareGroupsVCFOutputFormat {
    @BeforeClass
    public static void init() {

        GobyRengine.getInstance();

    }

    DiscoverVariantIterateSortedAlignments iterator;
    CompareGroupsVCFOutputFormat format;
    private int groupIndexA = 0;
    private int groupIndexB = 1;
    private int refIndex = 0;
    private int position = 0;
    private PrintWriter output;
    private DiscoverSequenceVariantsMode mode;
    private DifferentialExpressionAnalysis diffExpAnalysis;
    private DifferentialExpressionCalculator diffExpCalculator;
    private String[] groups = {"A", "B"};
    private String[] samples = {"s1", "s2", "s3", "s4"};
    private int[] readerIndexToGroupIndex;
    private VCFWriter statWriter;
    private int fisherExactPValueColumnIndex;

    @Before
    public void setUp() {

        iterator = createStrictMock(DiscoverVariantIterateSortedAlignments.class);
        mode = createMock(DiscoverSequenceVariantsMode.class);
        diffExpAnalysis = createMock(DifferentialExpressionAnalysis.class);
        diffExpCalculator = createMock(DifferentialExpressionCalculator.class);
        final StringWriter stringWriter = new StringWriter();
        output = new PrintWriter(stringWriter);
        format = new CompareGroupsVCFOutputFormat();
        statWriter = createMock("statWriter", VCFWriter.class);

        statWriter.clearAlternateAlleles();
        statWriter.setInfo(eq(2), anyDouble());
        statWriter.setInfo(eq(3), anyDouble());
        statWriter.setInfo(eq(4), anyDouble());


        expectLastCall().anyTimes();
        statWriter.setInfo(anyInt(), anyInt());
        expectLastCall().anyTimes();
        statWriter.setInfo(anyInt(), EasyMock.<CharSequence>anyObject());
        expectLastCall().anyTimes();
        statWriter.setSampleValue(anyInt(), anyInt(), EasyMock.<CharSequence>anyObject());
        expectLastCall().anyTimes();
        statWriter.setSampleValue(anyInt(), anyInt(), anyInt());

        expectLastCall().anyTimes();
        statWriter.setChromosome(EasyMock.<CharSequence>anyObject());
        statWriter.setPosition(position + 1); // position is written 1-based
        statWriter.setReferenceAllele(EasyMock.anyObject(String.class));
        expectLastCall().atLeastOnce();


        statWriter.addAlternateAllele(EasyMock.<String>anyObject());
        expectLastCall().atLeastOnce();

        expect(statWriter.codeGenotype(EasyMock.<String>anyObject())).andReturn(new MutableString("0/1/2"));
        expectLastCall().anyTimes();
        statWriter.setFlag(EasyMock.anyInt(), EasyMock.anyBoolean());
        expectLastCall().anyTimes();
        statWriter.writeRecord();
        expectLastCall().anyTimes();

        expect(mode.getDiffExpAnalyzer()).andReturn(diffExpAnalysis);
        expect(mode.getDiffExpCalculator()).andReturn(diffExpCalculator);
        expect(mode.getGroups()).andReturn(groups);
        expect(mode.getSamples()).andReturn(samples);
        expect(mode.getReadIndexStats()).andReturn(null);
    }


    @Test
    public void testNoDifference() throws Exception {
        synchronized (GobyRengine.getInstance().getRengine()) {
            SampleCountInfo[] sampleCounts = makeSampleCounts(1, 30, 5, 30, 6, 2, 60, 10, 60, 12);
            ObjectArrayList<PositionBaseInfo> list = new ObjectArrayList<PositionBaseInfo>();
            expect(mode.getReaderIndexToGroupIndex()).andReturn(readerIndexToGroupIndex);
            replay(mode);

            statWriter.setInfo(eq(fisherExactPValueColumnIndex = 5), eq(1d));

            replay(statWriter);
            format.allocateStorage(20, 2);
            format.defineColumns(output, mode);
            format.setStatWriter(statWriter);
            format.writeRecord(iterator, sampleCounts,
                    refIndex,
                    position,
                    list,
                    groupIndexA,
                    groupIndexB);
            verify(mode);
            verify(statWriter);
        }
    }


    @Test
    public void testAllelicDifference() throws Exception {
        synchronized (GobyRengine.getInstance().getRengine()) {
            SampleCountInfo[] sampleCounts = makeSampleCounts(1, 30, 5, 30, 6, 2 /* this value would have been 60 is no difference existed for this allele */, 0, 10, 60, 12);
            ObjectArrayList<PositionBaseInfo> list = new ObjectArrayList<PositionBaseInfo>();
            expect(mode.getReaderIndexToGroupIndex()).andReturn(readerIndexToGroupIndex);
            replay(mode);

            statWriter.setInfo(eq(fisherExactPValueColumnIndex = 5), lt(1d));
            replay(statWriter);
            format.allocateStorage(20, 2);
            format.defineColumns(output, mode);
            format.setStatWriter(statWriter);
            format.writeRecord(iterator, sampleCounts,
                    refIndex,
                    position,
                    list,
                    groupIndexA,
                    groupIndexB);
            verify(mode);
            verify(statWriter);
        }
    }

    @Test
    public void testAlignIndels() {
        SampleCountInfo[] sampleCounts = makeSampleCountsWithIndels();


        assertEquals(SampleCountInfo.BASE_MAX_INDEX + 2, sampleCounts[0].getGenotypeMaxIndex());
        assertEquals(SampleCountInfo.BASE_MAX_INDEX + 1, sampleCounts[1].getGenotypeMaxIndex());
        assertTrue(sampleCounts[0].hasIndels());
        assertTrue(sampleCounts[1].hasIndels());

        SampleCountInfo.alignIndels(sampleCounts);
        assertEquals(SampleCountInfo.BASE_MAX_INDEX + 2, sampleCounts[0].getGenotypeMaxIndex());
        assertEquals(SampleCountInfo.BASE_MAX_INDEX + 2, sampleCounts[1].getGenotypeMaxIndex());
        String sample_0_indel0_to = getGenotypeAsString(sampleCounts[0], SampleCountInfo.BASE_MAX_INDEX);
        String sample_0_indel1_to = getGenotypeAsString(sampleCounts[0], SampleCountInfo.BASE_MAX_INDEX + 1);
        String sample_1_indel0_to = getGenotypeAsString(sampleCounts[1], SampleCountInfo.BASE_MAX_INDEX);
        String sample_1_indel1_to = getGenotypeAsString(sampleCounts[1], SampleCountInfo.BASE_MAX_INDEX + 1);

        assertEquals("first indel genotype must match across both samples", sample_0_indel0_to, sample_1_indel0_to);
        assertEquals("second indel genotype must match across both samples", sample_0_indel1_to, sample_1_indel1_to);
        assertEquals(5, sampleCounts[0].getGenotypeCount(SampleCountInfo.BASE_MAX_INDEX));
        assertEquals(10, sampleCounts[0].getGenotypeCount(SampleCountInfo.BASE_MAX_INDEX + 1));
        assertEquals("frequency of dummy indel must be zero, dummy indels were never observed in the sample.",
                0, sampleCounts[1].getGenotypeCount(SampleCountInfo.BASE_MAX_INDEX));
        assertEquals(3, sampleCounts[1].getGenotypeCount(SampleCountInfo.BASE_MAX_INDEX + 1));
    }

    private String getGenotypeAsString(SampleCountInfo sampleCount, int genotypeIndex) {
        MutableString dest = new MutableString();
        sampleCount.getGenotype(genotypeIndex, dest);
        return dest.toString();
    }

    private SampleCountInfo[] makeSampleCountsWithIndels() {
        SampleCountInfo[] sampleCounts = makeSampleCounts(1, 30, 5, 30, 6, 2 /* this value would have been 60 is no difference existed for this allele */, 0, 10, 60, 12);

        // add indels to these sample counts, TODO:
        EquivalentIndelRegion indel1 = new EquivalentIndelRegion();
        indel1.startPosition = 1;
        indel1.endPosition = 5;
        indel1.referenceIndex = 0;
        indel1.from = "AC";
        indel1.to = "A---C";
        indel1.flankLeft = "T";
        indel1.flankRight = "";
        indel1.sampleIndex = 0;
        indel1.setFrequency(5);
        EquivalentIndelRegion indel2 = new EquivalentIndelRegion();
        indel2.startPosition = 1;
        indel2.endPosition = 3;
        indel2.referenceIndex = 0;
        indel2.from = "AC";
        indel2.to = "ATC";
        indel2.flankLeft = "T";
        indel2.flankRight = "";
        indel1.sampleIndex = 0;
        indel2.setFrequency(10);
        sampleCounts[0].addIndel(indel1);
        sampleCounts[0].addIndel(indel2);
        EquivalentIndelRegion indel3 = indel2.copy();
        indel3.sampleIndex = 1;
        indel3.setFrequency(3);

        sampleCounts[1].addIndel(indel3);
        return sampleCounts;
    }

    @Test
    public void testIndelDifferences() throws Exception {
        if (Release1_9_7_2.callIndels) {
            synchronized (GobyRengine.getInstance().getRengine()) {
                SampleCountInfo[] sampleCounts = makeSampleCountsWithIndels();

                ObjectArrayList<PositionBaseInfo> list = new ObjectArrayList<PositionBaseInfo>();
                expect(mode.getReaderIndexToGroupIndex()).andReturn(readerIndexToGroupIndex);
                replay(mode);

                statWriter.setInfo(eq(fisherExactPValueColumnIndex = 5), lt(1d));
                replay(statWriter);
                format.allocateStorage(20, 2);
                format.defineColumns(output, mode);
                format.setStatWriter(statWriter);
                SampleCountInfo.alignIndels(sampleCounts);
                format.writeRecord(iterator, sampleCounts,
                        refIndex,
                        position,
                        list,
                        groupIndexA,
                        groupIndexB);
                verify(mode);
                verify(statWriter);
            }
        }
    }

    private SampleCountInfo[] makeSampleCounts(int a_0, int t_0, int c_0, int refCount_0, int varCount_0,
                                               int a_1, int t_1, int c_1, int refCount_1, int varCount_1) {
        SampleCountInfo[] sampleCounts = new SampleCountInfo[20];
        int i = 0;
        readerIndexToGroupIndex = new int[20];
        for (i = 0; i < 10; i++) {
            sampleCounts[i] = new SampleCountInfo();
            sampleCounts[i].counts[SampleCountInfo.BASE_A_INDEX] = a_0;
            sampleCounts[i].counts[SampleCountInfo.BASE_T_INDEX] = t_0;
            sampleCounts[i].counts[SampleCountInfo.BASE_C_INDEX] = c_0;
            sampleCounts[i].counts[SampleCountInfo.BASE_OTHER_INDEX] = 0;
            sampleCounts[i].referenceBase = 'C';
            sampleCounts[i].refCount = refCount_0;
            sampleCounts[i].varCount = varCount_0;
            sampleCounts[i].sampleIndex = i;
            readerIndexToGroupIndex[i] = 0;
        }
        for (; i < 20; i++) {
            sampleCounts[i] = new SampleCountInfo();
            sampleCounts[i].counts[SampleCountInfo.BASE_A_INDEX] = a_1;
            sampleCounts[i].counts[SampleCountInfo.BASE_T_INDEX] = t_1;
            sampleCounts[i].counts[SampleCountInfo.BASE_C_INDEX] = c_1;
            sampleCounts[i].counts[SampleCountInfo.BASE_OTHER_INDEX] = 0;
            sampleCounts[i].referenceBase = 'C';
            sampleCounts[i].refCount = refCount_1;
            sampleCounts[i].varCount = varCount_1;
            sampleCounts[i].sampleIndex = i;
            readerIndexToGroupIndex[i] = 1;
        }

        return sampleCounts;
    }
}
