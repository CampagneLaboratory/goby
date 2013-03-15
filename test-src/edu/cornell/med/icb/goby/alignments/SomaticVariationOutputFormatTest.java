package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.alignments.processors.AlignmentProcessorInterface;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.junit.Test;
import org.easymock.EasyMock;

import static org.easymock.EasyMock.expect;
import static org.easymock.EasyMock.replay;

/**
 * @author Fabien Campagne
 *         Date: 3/13/13
 *         Time: 1:49 PM
 */
public class SomaticVariationOutputFormatTest {
    private static final int GERMLINE = 0;
    private static final int SOMATIC = 1;
    int[] dummyIndices = new int[]{-1, -1};

    @Test
    public void testEstimateSomaticPValue1() throws Exception {
        SampleCountInfo[] sampleCounts = new SampleCountInfo[2];
        sampleCounts[GERMLINE] = new SampleCountInfo();
        sampleCounts[SOMATIC] = new SampleCountInfo();

        assertPWithCounts(1.0, sampleCounts);

    }

    @Test
    public void testEstimateSomaticPValue2() throws Exception {
        SampleCountInfo[] sampleCounts = new SampleCountInfo[2];
        sampleCounts[GERMLINE] = new SampleCountInfo();
        sampleCounts[SOMATIC] = new SampleCountInfo();
        sampleCounts[GERMLINE].counts[SampleCountInfo.BASE_A_INDEX] = 0;
        sampleCounts[GERMLINE].counts[SampleCountInfo.BASE_C_INDEX] = 0;
        sampleCounts[SOMATIC].counts[SampleCountInfo.BASE_C_INDEX] = 100;
        sampleCounts[SOMATIC].counts[SampleCountInfo.BASE_A_INDEX] = 0;
        assertPWithCounts(2.2229697512078737E-19, sampleCounts);

    }

    @Test
    public void testEstimateSomaticPValue3() throws Exception {
        SampleCountInfo[] sampleCounts = new SampleCountInfo[2];
        sampleCounts[GERMLINE] = new SampleCountInfo();
        sampleCounts[SOMATIC] = new SampleCountInfo();
        sampleCounts[GERMLINE].counts[SampleCountInfo.BASE_A_INDEX] = 50;
        sampleCounts[GERMLINE].counts[SampleCountInfo.BASE_C_INDEX] = 0;
        sampleCounts[SOMATIC].counts[SampleCountInfo.BASE_C_INDEX] = 100;
        sampleCounts[SOMATIC].counts[SampleCountInfo.BASE_A_INDEX] = 50;
        assertPWithCounts(2.2229697512078737E-19, sampleCounts);

    }
    @Test
    public void testEstimateSomaticPValue4() throws Exception {
//        The following hit must not be found because they are more than 5 hits in the germline sample:
  //      1	28422629	.	TCCC	T-CC	.	.	BIOMART_COORDS=1:28422629-28422629;INDEL;Somatic-P-value(Fisher)[YHUDRHR-LMFRNGS-7-42-F-FSGSR-PBMC-patient-T0]=0.002235811068873165	GT:BC:GB:FB	0/1:T=388,T-CC=6:394:0	0/1:T=274,T-CC=6:280:0

        SampleCountInfo[] sampleCounts = new SampleCountInfo[2];
        sampleCounts[GERMLINE] = new SampleCountInfo();
        sampleCounts[SOMATIC] = new SampleCountInfo();
        sampleCounts[GERMLINE].counts[SampleCountInfo.BASE_A_INDEX] = 50;
        sampleCounts[GERMLINE].counts[SampleCountInfo.BASE_C_INDEX] = 5;     // too many bases for C genotype in germline. P-value will be 1.0.
        sampleCounts[SOMATIC].counts[SampleCountInfo.BASE_C_INDEX] = 100;
        sampleCounts[SOMATIC].counts[SampleCountInfo.BASE_A_INDEX] = 50;
        assertPWithCounts(1.0, sampleCounts);

    }

    private void assertPWithCounts(double expectedPValue, SampleCountInfo[] sampleCounts) {
        VCFWriter statsWriter = EasyMock.createMock(VCFWriter.class);
        sampleCounts[0].sampleIndex=0;
        sampleCounts[1].sampleIndex=1;
        SomaticVariationOutputFormat output = new SomaticVariationOutputFormat();
        output.setupR();
        output.allocateStorage(2,0);
        int[][] sampleIndex2GermlineSampleIndices = new int[][]{new int[]{-1}, new int[]{0}};
        IntArrayList indices = new IntArrayList();
        // sample 0 is germline, 1 is somatic:
        indices.add(1);
        output.setSomaticSampleIndices(indices);
        output.setSample2FatherSampleIndex(dummyIndices);
        output.setSample2MotherSampleIndex(dummyIndices);
        output.setSample2GermlineSampleIndices(sampleIndex2GermlineSampleIndices);


        statsWriter.setInfo(-1, expectedPValue);
        output.setSomaticPValueIndex(dummyIndices);
        output.setStatsWriter(statsWriter);
        replay(statsWriter);
        output.estimateSomaticPValue(sampleCounts);
    }
}
