/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.readers.sam;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;


/**
 * Test SamRecordParser.
 */
public class TestGobyPaperTop5000s {

    private static final Logger LOG = Logger.getLogger(TestGobyPaperTop5000s.class);

    private static final String BASE_TEST_INPUT_DIR = "test-data/goby-paper-top-5000s";
    private static final String BASE_TEST_OUTPUT_DIR = "test-results/goby-paper-top-5000s";

    /**
     * Setup the filenames for a RoundTripAlignment in a standard fashion for the 5000 tests.
     * @param input the input sam/bam filename. Should end with .sam, .bam, or .sam.gz
     * @return a RoundTripAlignment with input and output filenames configured
     */
    private RoundTripAlignment setupRoundTrip(final String input) {
        final String basename;
        if (input.toLowerCase().endsWith(".gz")) {
            basename = FilenameUtils.getBaseName(FilenameUtils.getBaseName(input));
        } else {
            basename = FilenameUtils.getBaseName(input);
        }
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.sourceBamFilename = FilenameUtils.concat(BASE_TEST_INPUT_DIR, input);
        rtc.destGobyBasename = FilenameUtils.concat(BASE_TEST_OUTPUT_DIR, basename);
        rtc.destBamFilename = FilenameUtils.concat(BASE_TEST_OUTPUT_DIR, basename + ".sam");
        return rtc;
    }

    /**
     * The first 5000 alignments of HZFWPTI. Round trip test.
     *
     * @throws java.io.IOException error
     */
    @Test
    public void testRoundTripFirst5000HZFWPTI() throws IOException {
        final RoundTripAlignment rtc = setupRoundTrip("HZFWPTI-5000.sam.gz");
        rtc.inputGenomeFilename = findThousandGenome();
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
        rtc.keepSoftClips = false;
        rtc.testRoundTripAny();
    }

    /**
     * The first 5000 alignments of UANMNXR. Round trip test.
     *
     * @throws java.io.IOException error
     */
    @Test
    public void testRoundTripFirst5000UANMNXR() throws IOException {
        final RoundTripAlignment rtc = setupRoundTrip("UANMNXR-5000.sam.gz");
        rtc.inputGenomeFilename = findThousandGenome();
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * The first 5000 alignments of MYHZZJH. Round trip test.
     * Spliced.
     *
     * NOTE: There are 11 reads where the source read has an N and the destination read as a non-N base.
     * @throws java.io.IOException error
     */
    @Test
    public void testRoundTripFirst5000MYHZZJH() throws IOException {
        final RoundTripAlignment rtc = setupRoundTrip("MYHZZJH-5000.sam.gz");
        rtc.inputGenomeFilename = findThousandGenome();
        rtc.sortGoby = true;
        rtc.sortSam = true;
        rtc.compareAtEachNewPosition = false;
        rtc.allowSourceNs = true;
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * The first 5000 alignments of ZHUUJKS. Round trip test.
     * Spliced.
     *
     * NOTE: There are 60 reads where the source read has an N and the destination read as a non-N base.
     * @throws java.io.IOException error
     */
    @Test
    public void testRoundTripFirst5000ZHUUJKS() throws IOException {
        final RoundTripAlignment rtc = setupRoundTrip("ZHUUJKS-5000.sam.gz");
        rtc.inputGenomeFilename = findThousandGenome();
        rtc.sortGoby = true;
        rtc.sortSam = true;
        rtc.compareAtEachNewPosition = false;
        rtc.allowSourceNs = true;
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * The first 5000 alignments of EJOYQAZ. Round trip test.
     * Spliced with tophat.
     * @throws java.io.IOException error
     */
    @Test
    public void testRoundTripFirst5000EJOYQAZ() throws IOException {
        final RoundTripAlignment rtc = setupRoundTrip("EJOYQAZ-5000.sam.gz");
        rtc.inputGenomeFilename = findHG19();
        rtc.sortGoby = true;
        rtc.sortSam = true;
        rtc.compareAtEachNewPosition = false;
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * The first 5000 alignments of JRO. Round trip test.
     * @throws java.io.IOException error
     */
    @Test
    public void testRoundTripFirst5000JRODTYG() throws IOException {
        final RoundTripAlignment rtc = setupRoundTrip("JRODTYG-5000.sam.gz");
        rtc.inputGenomeFilename = findMM9();
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * The first 5000 alignments of ZVLRRJH. Round trip test.
     * @throws java.io.IOException error
     */
    @Test
    public void testRoundTripFirst5000ZVLRRJH() throws IOException {
        final RoundTripAlignment rtc = setupRoundTrip("ZVLRRJH-5000.sam.gz");
        rtc.inputGenomeFilename = findHG18();
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * The first 5000 alignments of XAAOBVT. Round trip test.
     * @throws java.io.IOException error
     */
    @Test
    public void testRoundTripFirst5000XAAOBVT() throws IOException {
        final RoundTripAlignment rtc = setupRoundTrip("XAAOBVT-5000.sam.gz");
        rtc.inputGenomeFilename = findHG18();
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * The first 5000 alignments of UCCWRUX. Round trip test.
     * @throws java.io.IOException error
     */
    @Test
    public void testRoundTripFirst5000UCCWRUX() throws IOException {
        final RoundTripAlignment rtc = setupRoundTrip("UCCWRUX-5000.sam.gz");
        rtc.inputGenomeFilename = findHG18();
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * The first 5000 alignments of HENGLIT. Round trip test.
     * @throws java.io.IOException error
     */
    @Test
    public void testRoundTripFirst5000HENGLIT() throws IOException {
        final RoundTripAlignment rtc = setupRoundTrip("HENGLIT-5000.sam.gz");
        rtc.inputGenomeFilename = findCelegansGenome();
        rtc.testRoundTripAny();
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * Full check of HZF. This dataset is NOT on the server so this test shouldn't be run on the server.
     *
     * @throws java.io.IOException error
     */
    // @Test
    public void testHzFullCompare() throws IOException {
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.inputGenomeFilename = findThousandGenome();
        rtc.sourceBamFilename = "/tmp/HZ-bam-bam/HZFWPTI-source.bam";
        rtc.destGobyBasename = "/tmp/HZ-bam-bam/HZFWPTI";
        rtc.destBamFilename = "/tmp/HZ-bam-bam/HZFWPTI.bam";
        rtc.keepQualityScores = true;
        rtc.keepSoftClips = true;
        rtc.convertBamToGoby = false;
        rtc.convertGobyToBam = false;
        rtc.testRoundTripAny();
    }

    /**
     * NOT FOR SERVER
     *
     * @throws java.io.IOException error
     */
    // @Test
    public void testFirst5000_HENGLIT_SOURCE_CRAM2() throws IOException {
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.sourceBamFilename = "/tmp/HENGLIT-source-5000.sam";
        rtc.destGobyBasename = "/tmp/HENGLIT-hybrid";
        rtc.destBamFilename = "/tmp/HENGLIT-cram2-5000.sam";
        rtc.convertBamToGoby = false;
        rtc.convertGobyToBam = false;
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * NOT FOR SERVER
     *
     * @throws java.io.IOException error
     */
    // @Test
    public void testFirst5000_HZFWPTI_SOURCE_GOBY() throws IOException {
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.sourceBamFilename = "/tmp/HZFWPTI-source-5000.sam";
        rtc.destGobyBasename = "/tmp/HZFWPTI-hybrid";
        rtc.destBamFilename = "/tmp/HZFWPTI-top-5000-goby.sam";
        rtc.convertBamToGoby = false;
        rtc.convertGobyToBam = false;
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * NOT FOR SERVER
     *
     * @throws java.io.IOException error
     */
    // @Test
    public void testFirst5000_HZFWPTI_SOURCE_CRAM2() throws IOException {
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.sourceBamFilename = "/tmp/HZFWPTI-source-5000.sam";
        rtc.destGobyBasename = "/tmp/HZFWPTI-hybrid";
        rtc.destBamFilename = "/tmp/HZFWPTI-top-5000-cram2.sam";
        rtc.convertBamToGoby = false;
        rtc.convertGobyToBam = false;
        rtc.keepQualityScores = false;
        rtc.testRoundTripAny();
    }

    /**
     * The first 1M alignments of UAN. This large dataset does not exist on the testing server.
     *
     * @throws java.io.IOException error
     */
    // @Test
    public void testRoundTrip1M() throws IOException {
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.inputGenomeFilename = findThousandGenome();
        rtc.sourceBamFilename = "test-data/splicedsamhelper/1M.bam";
        rtc.destGobyBasename = FilenameUtils.concat(BASE_TEST_OUTPUT_DIR, "1M");
        rtc.destBamFilename = FilenameUtils.concat(BASE_TEST_OUTPUT_DIR, "1M.bam");
        rtc.canonicalMdzForComparison = false;
        final long start = System.currentTimeMillis();
        rtc.testRoundTripAny();
        System.out.printf("Execution time in ms=%d%n", System.currentTimeMillis() - start);
    }

    /**
     * The whole of HENGLIT. This large dataset does not exist on the testing server.
     *
     * @throws java.io.IOException error
     */
    // @Test
    public void testRoundTripHenglit() throws IOException {
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.inputGenomeFilename = findCelegansGenome();
        rtc.sourceBamFilename = "/tmp/HENGLIT.bam";
        rtc.destGobyBasename = "/tmp/HENGLIT-to-goby";
        rtc.destBamFilename = "/tmp/HENGLIT-from-goby.bam";
        //rtc.convertBamToGoby = false;
        //rtc.convertGobyToBam = false;
        rtc.canonicalMdzForComparison = false;
        rtc.testRoundTripAny();
    }

    /**
     * The whole of HENGLIT. This large dataset does not exist on the testing server.
     * THIS version doesn't store read quals to better test the comparison when read quals aren't stored.
     *
     * @throws java.io.IOException error
     */
    // @Test
    public void testRoundTripHenglitNoQuals() throws IOException {
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.inputGenomeFilename = findCelegansGenome();
        rtc.sourceBamFilename = "/tmp/HENGLIT.bam";
        rtc.destGobyBasename = "/tmp/HENGLIT-to-goby-no-qual";
        rtc.destBamFilename = "/tmp/HENGLIT-from-goby-no-qual.bam";
        //rtc.convertBamToGoby = false;
        //rtc.convertGobyToBam = false;
        rtc.keepQualityScores = false;
        rtc.canonicalMdzForComparison = false;
        rtc.testRoundTripAny();
    }

    /**
     * Test that SHOULD fail, for local testing, not for server testing.
     * This is designed to be run manually as the 1M doesn't exist on the testing server.
     * Before running this test, run testRoundTrip1M().
     *
     * @throws java.io.IOException error
     */
    // @Test
    public void testRoundTripFail() throws IOException {
        final RoundTripAlignment rtc = new RoundTripAlignment();
        rtc.sourceBamFilename = FilenameUtils.concat(BASE_TEST_INPUT_DIR, "HZFWPTI-5000.sam.gz");
        rtc.destBamFilename = FilenameUtils.concat(BASE_TEST_INPUT_DIR, "JRODTYG-5000.sam.gz");
        rtc.convertBamToGoby = false;
        rtc.convertGobyToBam = false;
        rtc.testRoundTripAny();
    }

    private static String findGenome(final String[] dirs) {
        for (final String dir : dirs) {
            final String testRootFilename = dir + "/" + "random-access-genome";
            final String testFilename = testRootFilename + ".names";
            System.out.println("Looking for :" + testFilename);
            final File testFile = new File(testFilename);
            if (testFile.exists()) {
                return testRootFilename;
            }
        }
        return null;
    }

    /**
     * Find the local celegans random genome from an array of options.
     *
     * @return the local celegans random genome.
     */
    public static String findCelegansGenome() {
        // The following directories will be scanned in sequence to find a locatio where the genome is installed. Add new
        // locations as needed to indicate where the genome is installed on your local machine.
        final String[] dirs = {
                "/tmp/celegans",
                "/scratchLocal/gobyweb/input-data/reference-db/goby-benchmark-paper/cElegans",
                "/home/ccontrol/goby-data/celegans",
                "/data/cElegans"
        };
        return findGenome(dirs);
    }


    /**
     * Find the local 1000g random genome from an array of options.
     *
     * @return the local random genome.
     */
    public static String findThousandGenome() {
        final String[] dirs = {
                "/tmp/1000g",
                "/scratchLocal/gobyweb/input-data/reference-db/1000GENOMES.37/homo_sapiens/reference",
                "/home/ccontrol/goby-data/1000g",
                "/data/1000g"};
        return findGenome(dirs);
    }

    /**
     * Find the local HG19 random genome from an array of options.
     *
     * @return the random genome.
     */
    public static String findHG19() {
        final String[] dirs = {
                "/tmp/hg19",
                "/scratchLocal/gobyweb/input-data/reference-db/goby-benchmark-paper/hg19",
                "/home/ccontrol/goby-data/hg19",
                "/data/hg19"};
        return findGenome(dirs);
    }

    /**
     * Find the local HG19 random genome from an array of options.
     *
     * @return the random genome.
     */
    public static String findHG18() {
        final String[] dirs = {
                "/tmp/hg18",
                "/scratchLocal/gobyweb/input-data/reference-db/goby-benchmark-paper/hg18",
                "/home/ccontrol/goby-data/hg18",
                "/data/hg18"};
        return findGenome(dirs);
    }

    /**
     * Find the local mm9 random genome from an array of options.
     *
     * @return the local mm9  random genome.
     */
    public static String findMM9() {
        final String[] dirs = {
                "/tmp/mm9",
                "/scratchLocal/gobyweb/input-data/reference-db/goby-benchmark-paper/mm9",
                "/home/ccontrol/goby-data/mm9",
                "/data/mm9"
        };
        return findGenome(dirs);
    }

    @BeforeClass
    public static void beforeClass() {
        new File(BASE_TEST_OUTPUT_DIR).mkdirs();
    }
}
