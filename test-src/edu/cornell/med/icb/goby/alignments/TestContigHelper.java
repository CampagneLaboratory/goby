package edu.cornell.med.icb.goby.alignments;

import org.junit.Test;
import static junit.framework.Assert.assertEquals;
import java.io.File;

/**
 *
 * User: Eric Minwei Liu
 * Date: 11/26/13
 * Time: 3:57 PM
 *.
 */
public class TestContigHelper {

    // test case
    // 0__len__721:203-203	-	chr2	9630375  9630375
    // 1__len__681:107-107	+	chr10	92672164  92672164

    //aaaaa

    String contigName1 = "0__len__721";
    int contigPos1 = 203;

    String contigName2 = "1__len__681";
    int contigPos2 = 107;

    // non-exist contig id in the alntable
    String contigName3 = "0__len__72111";
    int contigPos3 = 204;



    @Test
    public void testTranslatedPosCase1() throws Exception {
        ContigHelper contigMapper = new ContigHelper(new File("test-data/contig-mapping/aln_mix_to_hg19.pos_mapping"));
        String testContigChr1 = contigMapper.translateContigID(contigName1);
        int testContigPos1 = contigMapper.translateContigPosition(contigName1, contigPos1);

        assertEquals("contigName should be equal", "chr2", testContigChr1 );
        assertEquals("contigPos should be equal", 9630375, testContigPos1 );

    }

    @Test
    public void testTranslatedPosCase2() throws Exception {
        ContigHelper contigMapper = new ContigHelper(new File("test-data/contig-mapping/aln_mix_to_hg19.pos_mapping"));
        String testContigChr2 = contigMapper.translateContigID(contigName2);
        int testContigPos2 = contigMapper.translateContigPosition(contigName2, contigPos2);

        assertEquals("contigName should be equal", "chr10", testContigChr2 );
        assertEquals("contigPos should be equal", 92672164, testContigPos2 );

    }

    @Test
    public void testTranslatedPosCase3() throws Exception {
        ContigHelper contigMapper = new ContigHelper(new File("test-data/contig_mapping/aln_mix_to_hg19.pos_mapping"));
        String testContigChr3 = contigMapper.translateContigID(contigName3);
        int testContigPos3 = contigMapper.translateContigPosition(contigName3, contigPos3);

        assertEquals("contigName should be equal", "0__len__72111", testContigChr3 );
        assertEquals("contigPos should be equal", 204, testContigPos3 );

    }

}
