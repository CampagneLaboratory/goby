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

import edu.cornell.med.icb.goby.alignments.DiscoverVariantIterateSortedAlignments;
import edu.cornell.med.icb.goby.alignments.PositionBaseInfo;
import edu.cornell.med.icb.goby.alignments.SampleCountInfo;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.junit.Before;
import org.junit.Test;

import java.io.PrintWriter;
import java.io.StringWriter;

import static org.easymock.EasyMock.*;

/**
 * @author Fabien Campagne
 *         Date: 5/30/11
 *         Time: 10:37 AM
 */
@SuppressWarnings({"IOResourceOpenedButNotSafelyClosed"})
public class TestCompareGroupsVCFOutputFormat {
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
    private String[] groups= {"A","B"};
    private String[] samples= {"s1","s2","s3","s4"};
    private int[] readerIndexToGroupIndex={0,0,1,1};
    @Before
    public void setUp() {

        iterator = createStrictMock(DiscoverVariantIterateSortedAlignments.class);
        mode = createMock(DiscoverSequenceVariantsMode.class);
        diffExpAnalysis=createMock(DifferentialExpressionAnalysis.class);
        diffExpCalculator=createMock(DifferentialExpressionCalculator.class);
        final StringWriter stringWriter = new StringWriter();
        output = new PrintWriter(stringWriter);
        format = new CompareGroupsVCFOutputFormat();


    }

    @Test
    public void testWriteRecord() throws Exception {
        SampleCountInfo[] sampleCounts = new SampleCountInfo[2];
        sampleCounts[0]=new SampleCountInfo();
        sampleCounts[1]=new SampleCountInfo();
        ObjectArrayList<PositionBaseInfo> list = new ObjectArrayList<PositionBaseInfo>();
        expect(mode.getDiffExpAnalyzer()).andReturn(diffExpAnalysis);
        expect(mode.getDiffExpCalculator()).andReturn(diffExpCalculator);
        expect(mode.getGroups()).andReturn(groups);
        expect(mode.getSamples()).andReturn(samples);
        expect(mode.getReaderIndexToGroupIndex()).andReturn(readerIndexToGroupIndex);
        expect(mode.getReadIndexStats()).andReturn(null);
        replay(mode);
        format.allocateStorage(20, 2);
        format.defineColumns(output, mode);
        format.writeRecord(iterator, sampleCounts,
                refIndex,
                position,
                list,
                groupIndexA,
                groupIndexB);
        verify(mode);
    }
}
