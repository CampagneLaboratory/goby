package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 3/15/13
 *         Time: 3:58 PM
 */
public class CommonIndelArtifactFilterTest {
    @Test
    public void testCountRepetitiveBases() throws Exception {
        CommonIndelArtifactFilter filter = new CommonIndelArtifactFilter();
        assertEquals(3, filter.countRepetitiveBases("AAAA"));
        assertEquals(0, filter.countRepetitiveBases("A"));
        assertEquals(3, filter.countRepetitiveBases("TAAAA"));
        assertEquals(3, filter.countRepetitiveBases("TAAAAG"));
        assertEquals(3, filter.countRepetitiveBases("TAAAGAAG"));

        assertEquals(6, filter.countRepetitiveBases("TUGGGGGGG"));
        assertEquals(0, filter.countRepetitiveBases(""));
        assertEquals(0, filter.countRepetitiveBases("-----A"));
        assertEquals(1, filter.countRepetitiveBases("-ATCTGCAA"));


    }
    @Test
       public void testCase1() throws Exception {
           CommonIndelArtifactFilter filter = new CommonIndelArtifactFilter();

        EquivalentIndelRegion indel = new EquivalentIndelRegion();
        indel.from="C-TTTTTTTTTTTTTT";
        indel.to="CTTTTTTTTTTTTTTT";
        assertEquals(14, filter.countRepetitiveBases(indel));



       }

}
