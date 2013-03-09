package edu.cornell.med.icb.goby.algorithmic.data;

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 3/9/13
 *         Time: 11:38 AM
 */
public class CovariateInfoTest {
    @Test
    public void testLoad() throws IOException {
        CovariateInfo info = CovariateInfo.parse("test-data/covariates/example-1.tsv");
        ObjectArraySet<String> expectedMales = new ObjectArraySet<String>();
        expectedMales.add("S1");
        expectedMales.add("S3");
        expectedMales.add("S4");
        assertEquals(expectedMales, info.samplesWithExactCovariate("gender", "Male"));

        ObjectArraySet<String> expectedSomatic = new ObjectArraySet<String>();
        expectedSomatic.add("S3");
        assertEquals(expectedSomatic, info.samplesWithExactCovariate("kind-of-sample", "Somatic"));
    }

    @Test
    public void testQueryContains() throws IOException {
        CovariateInfo info = CovariateInfo.parse("test-data/covariates/example-2.tsv");
        ObjectArraySet<String> expectedParentsOfS1 = new ObjectArraySet<String>();

        expectedParentsOfS1.add("S3");
        expectedParentsOfS1.add("S4");
        assertEquals(expectedParentsOfS1, info.samplesContainCovariate("parents", "S1"));

        assertEquals("S1|S2", info.getCovariateValue("S3","parents"));


        assertEquals("sample does not exist", null, info.getCovariateValue("S--","parents"));

        assertEquals("covariate does not exist", null, info.getCovariateValue("S3","parent---s"));

    }

}
