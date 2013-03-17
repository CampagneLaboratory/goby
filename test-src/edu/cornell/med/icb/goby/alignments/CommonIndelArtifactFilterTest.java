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
        indel.from = "C-TTTTTTTTTTTTTT";
        indel.to = "CTTTTTTTTTTTTTTT";
        assertEquals(14, filter.countRepetitiveBases(indel));


    }

    @Test
    public void testGapLength() throws Exception {
        CommonIndelArtifactFilter filter = new CommonIndelArtifactFilter();

        assertEquals(4, filter.gapLength("A----G"));
        assertEquals(4, filter.gapLength("-A---G"));
        assertEquals(4, filter.gapLength("A---G-"));
        assertEquals(4, filter.gapLength("----"));
        assertEquals(4, filter.gapLength("A-C-T-G-A"));
    }

    @Test
    public void testRepetitivePatterns() throws Exception {
        CommonIndelArtifactFilter filter = new CommonIndelArtifactFilter();

        assertEquals(2, filter.repeatPatternLength("--AGAGAG", "AGAGAGAG"));
        assertEquals(3, filter.repeatPatternLength("---AGTAGTAGT", "AGTAGTAGT"));
        assertEquals(4, filter.repeatPatternLength("----CGATCGATCGAT ", "CGATCGAT"));


        EquivalentIndelRegion indel = new EquivalentIndelRegion();
        indel.from = "--AGAGAG";
        indel.to = "AGAGAGAG";
        assertEquals(2, filter.repeatPatternLength(indel));

        indel.from = "AGAGAGAG";
        indel.to = "--AGAGAG";
        assertEquals(2, filter.repeatPatternLength(indel));
    }
}
