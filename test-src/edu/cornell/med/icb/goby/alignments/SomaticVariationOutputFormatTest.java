package edu.cornell.med.icb.goby.alignments;

import com.martiansoftware.jsap.JSAPException;
import edu.cornell.med.icb.goby.algorithmic.data.CovariateInfo;
import edu.cornell.med.icb.goby.alignments.processors.AlignmentProcessorInterface;
import edu.cornell.med.icb.goby.stats.VCFWriter;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.junit.Test;
import org.easymock.EasyMock;

import java.io.IOException;

import static junit.framework.Assert.assertTrue;
import static org.easymock.EasyMock.expect;
import static org.easymock.EasyMock.replay;

/**
 * EasyMock for somatic variation calls.
 *
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

        assertPWithCounts(1.0, sampleCounts, 0f);

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
        assertPWithCounts(2.2229697512078737E-19, sampleCounts,100f);

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
        assertPWithCounts(2.2229697512078737E-19, sampleCounts, 66.66667f);

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
        assertPWithCounts(1.0, sampleCounts, 0f);

    }

    private void assertPWithCounts(double expectedPValue, SampleCountInfo[] sampleCounts, float expectedSomaticFrequency) {
        VCFWriter statsWriter = EasyMock.createMock(VCFWriter.class);
        sampleCounts[0].sampleIndex = 0;
        sampleCounts[1].sampleIndex = 1;
        SomaticVariationOutputFormat output = new SomaticVariationOutputFormat();
        output.setupR();
        output.allocateStorage(2, 0);
        int[][] sampleIndex2GermlineSampleIndices = new int[][]{new int[]{-1}, new int[]{0}};
        IntArrayList indices = new IntArrayList();
        // sample 0 is germline, 1 is somatic:
        indices.add(1);
        output.setSomaticSampleIndices(indices);
        output.setSample2FatherSampleIndex(dummyIndices);
        output.setSample2MotherSampleIndex(dummyIndices);
        output.setSample2GermlineSampleIndices(sampleIndex2GermlineSampleIndices);
        output.setCandidateFrequencyIndex(new int[]{1, 2, 3});

        statsWriter.setInfo(-1, expectedPValue);
        statsWriter.setInfo(2, expectedSomaticFrequency);
        output.setSomaticPValueIndex(dummyIndices);
        output.setStatsWriter(statsWriter);
        replay(statsWriter);
        output.allocateIsSomaticCandidate(sampleCounts);
        output.estimateSomaticPValue(sampleCounts);
    }


    @Test

    public void testDiscoverSomaticVariation2() throws IOException, JSAPException {
        SomaticVariationOutputFormat output = new SomaticVariationOutputFormat();
        output.setupR();
        output.allocateStorage(3, 0);

        final CovariateInfo covInfo = CovariateInfo.parse("test-data/covariates/example-4.tsv");
        output.setCovariateInfo(covInfo);
        IntArrayList indices = new IntArrayList();

        indices.add(1);
        VCFWriter statsWriter = EasyMock.createMock(VCFWriter.class);
        output.setStatsWriter(statsWriter);
        output.setSomaticSampleIndices(indices);
        // sample 0 is mother, 1 is somatic, 2 is father:
        output.setSample2FatherSampleIndex(new int[]{-1, 2, -1});
        output.setSample2MotherSampleIndex(new int[]{-1, 0, -1});
        // no germline in this trio test case, only somatic, father and mother.
        int[][] sampleIndex2GermlineSampleIndices = new int[][]{new int[]{-1}, new int[]{-1}, new int[]{-1},};
        output.setSample2GermlineSampleIndices(sampleIndex2GermlineSampleIndices);
        output.setCandidateFrequencyIndex(new int[]{2, 3, 4});
        output.setSomaticPValueIndex(new int[]{1, 2, 3});


        int somaticCount = 10;
        int somaticFailedCount = 0;
        int fatherCount = 0;
        int motherCount = 0;
        int baseline = 50;

        SampleCountInfo[] sampleCounts = makeSomaticSampleCounts(somaticCount, somaticFailedCount,
                fatherCount, motherCount, baseline);
        output.allocateIsSomaticCandidate(sampleCounts);


        assertTrue(output.isPossibleSomaticVariation(sampleCounts));
        statsWriter.setInfo(2, 0.016253869969040255d);
        statsWriter.setInfo(3, 9.090909f);
        replay(statsWriter);
        output.estimateSomaticPValue(sampleCounts);
        somaticCount = 20;
        sampleCounts=makeSomaticSampleCounts(somaticCount, somaticFailedCount,
                fatherCount, motherCount, baseline);
        output.allocateIsSomaticCandidate(sampleCounts);
        statsWriter = EasyMock.createMock(VCFWriter.class);
        output.setStatsWriter(statsWriter);
        statsWriter.setInfo(2, 2.1795989537925028E-4);
        statsWriter.setInfo(3, 16.666668f);
        replay(statsWriter);
        output.estimateSomaticPValue(sampleCounts);
    }

    private SampleCountInfo[] makeSomaticSampleCounts(int somaticCount, int somaticFailedCount,
                                                      int fatherCount, int motherCount, int baseline) {
        SampleCountInfo[] sampleCounts = new SampleCountInfo[3];
        sampleCounts[0] = new SampleCountInfo();
        sampleCounts[0].counts[SampleCountInfo.BASE_A_INDEX] = baseline;
        sampleCounts[0].counts[SampleCountInfo.BASE_T_INDEX] = motherCount;
        sampleCounts[0].counts[SampleCountInfo.BASE_C_INDEX] = baseline;
        sampleCounts[0].counts[SampleCountInfo.BASE_OTHER_INDEX] = 0;
        sampleCounts[0].referenceBase = 'A';
        sampleCounts[0].refCount = baseline * 2;
        sampleCounts[0].varCount = somaticCount;
        sampleCounts[0].sampleIndex = 0;
        sampleCounts[0].failedCount = somaticFailedCount;


        sampleCounts[1] = new SampleCountInfo();
        sampleCounts[1].counts[SampleCountInfo.BASE_A_INDEX] = baseline;
        sampleCounts[1].counts[SampleCountInfo.BASE_T_INDEX] = somaticCount;
        sampleCounts[1].counts[SampleCountInfo.BASE_C_INDEX] = baseline;
        sampleCounts[1].counts[SampleCountInfo.BASE_OTHER_INDEX] = 0;
        sampleCounts[1].referenceBase = 'A';
        sampleCounts[1].refCount = baseline * 2 + motherCount;
        sampleCounts[1].varCount = motherCount;
        sampleCounts[1].failedCount = somaticFailedCount;
        sampleCounts[1].sampleIndex = 1;

        sampleCounts[2] = new SampleCountInfo();
        sampleCounts[2].counts[SampleCountInfo.BASE_A_INDEX] = baseline;
        sampleCounts[2].counts[SampleCountInfo.BASE_T_INDEX] = fatherCount;
        sampleCounts[2].counts[SampleCountInfo.BASE_C_INDEX] = baseline;
        sampleCounts[2].counts[SampleCountInfo.BASE_OTHER_INDEX] = 0;
        sampleCounts[2].referenceBase = 'A';
        sampleCounts[2].refCount = baseline * 2 + fatherCount;
        sampleCounts[2].varCount = fatherCount;
        sampleCounts[2].failedCount = somaticFailedCount;
        sampleCounts[2].sampleIndex = 2;
        return sampleCounts;
    }

}
