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

package edu.cornell.med.icb.goby.stats;

import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.util.TestFiles;
import static junit.framework.Assert.assertTrue;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Test that StatsWriter can write both TSV and VCF files.
 *
 * @author Fabien Campagne
 *         Date: Mar 30, 2011
 *         Time: 6:58:10 PM
 */
public class TestStatsWriter extends TestFiles {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestStatsWriter.class);

    private static final String BASE_TEST_DIR = "test-results/stats-writer";

    @Test
    public void testTSV() throws IOException {

        final String file = FilenameUtils.concat(BASE_TEST_DIR, "1.tsv");

        TSVWriter writer = new TSVWriter(new PrintWriter(new FileWriter(new File(file))));

        writer.defineColumn("c1");
        writer.defineColumnAttributes("C", 1, ColumnType.String, "something about c1", "c1");
        writer.defineColumn("C2");
        writer.defineColumnAttributes("D", 1, ColumnType.Flag, "something about c2", "C2");
        writer.writeHeader();
        writer.close();

        assertEquals(new File("test-data/stats-writer/1.tsv"), new File("test-results/stats-writer/1.tsv"));
    }

    @Test
    public void testVCF() throws IOException {

        final String file = FilenameUtils.concat(BASE_TEST_DIR, "1.vcf");

        VCFWriter writer = new VCFWriter(new PrintWriter(new FileWriter(file)));

        int fieldC = writer.defineField("INFO", "C", 1, ColumnType.String, "something about C");
        int fieldD = writer.defineField("INFO", "D", 1, ColumnType.String, "something about D");

        writer.writeHeader();
        writer.setInfo(fieldC, "1");
        writer.setInfo(fieldD, "dvalue");
        writer.writeRecord();
        writer.close();

        assertEquals(new File("test-data/stats-writer/1.vcf"), new File("test-results/stats-writer/1.vcf"));
    }


    @Test
    public void testVCFFormat() throws IOException {

        final String file = FilenameUtils.concat(BASE_TEST_DIR, "format-1.vcf");

        VCFWriter writer = new VCFWriter(new PrintWriter(new FileWriter(file)));

        int fieldC = writer.defineField("FORMAT", "GT", 1, ColumnType.String, "Desc GT");
        int fieldD = writer.defineField("FORMAT", "GQ", 1, ColumnType.String, "Desc GQ");
        String samples[] = {"sample-id-1", "sample-id-2"};

        writer.defineSamples(samples);
        writer.writeHeader();
        writer.setSampleValue(fieldC, 0, "A");
        writer.setSampleValue(fieldD, 0, "B");
        writer.setSampleValue(fieldC, 1, "B");
        writer.setSampleValue(fieldD, 1, "A");

        writer.writeRecord();
        writer.close();

        assertEquals(new File("test-data/stats-writer/expected-format-1.vcf"), new File("test-results/stats-writer/format-1.vcf"));
    }
    @Test
        public void testVCFFormat2() throws IOException {

            final String file = FilenameUtils.concat(BASE_TEST_DIR, "format-2.vcf");

            VCFWriter writer = new VCFWriter(new PrintWriter(new FileWriter(file)));

            int fieldC = writer.defineField("FORMAT", "GT", 1, ColumnType.String, "Desc GT");
            int fieldD = writer.defineField("FORMAT", "GQ", 1, ColumnType.String, "Desc GQ");
            String samples[] = {"sample-id-1", "sample-id-2"};

            writer.defineSamples(samples);
            writer.writeHeader();
            writer.setSampleValue(fieldC, 0, "A");
            writer.setSampleValue(fieldD, 0, "B");

            writer.setSampleValue(fieldC, 1, "A");
               //   writer.setSampleValue(fieldD, 1, "B");

            writer.writeRecord();
            writer.close();

            assertEquals(new File("test-data/stats-writer/expected-format-2.vcf"), new File("test-results/stats-writer/format-2.vcf"));
        }

    @Test
    public void testCodeGenotype() throws IOException {
        final String file = FilenameUtils.concat(BASE_TEST_DIR, "not-tested.vcf");

        VCFWriter writer = new VCFWriter(new PrintWriter(new FileWriter(file)));
        writer.setReferenceAllele("A");
        // Only CC ref should remain in output
        writer.setReferenceAllele("CC");
        writer.addAlternateAllele("T");
        boolean success;
        try {
            writer.codeGenotype("A/CC/T");
            success = false;
        } catch (IllegalArgumentException e) {
            success = true;
        }
        assertTrue("incorrect alleles must raise IllegalArgumentException", success);
    }

    @Test
    public void testGenotypes() throws IOException {

        final String file = FilenameUtils.concat(BASE_TEST_DIR, "genotypes-1.vcf");

        VCFWriter writer = new VCFWriter(new PrintWriter(new FileWriter(file)));

        int fieldC = writer.defineField("FORMAT", "GT", 1, ColumnType.String, "Genotype");
        String samples[] = {"sample-id-1", "sample-id-2"};

        writer.defineSamples(samples);
        writer.writeHeader();
        writer.setReferenceAllele("A");
        // Only CC ref should remain in output
        writer.setReferenceAllele("CC");

        writer.addAlternateAllele("A");
        writer.addAlternateAllele("N");
        writer.addAlternateAllele("T");
        writer.addAlternateAllele("C");
        writer.setSampleValue(fieldC, 0, writer.codeGenotype("CC/T"));
        writer.setSampleValue(fieldC, 1, writer.codeGenotype("A/C/T"));

        writer.writeRecord();
        writer.setReferenceAllele("A");
        writer.addAlternateAllele("C");
        writer.addAlternateAllele("T");
        writer.setSampleValue(fieldC, 0, writer.codeGenotype("A/A"));
        writer.setSampleValue(fieldC, 1, writer.codeGenotype("A/C"));


        writer.writeRecord();
        writer.close();

        assertEquals(new File("test-data/stats-writer/expected-genotypes-1.vcf"), new File("test-results/stats-writer/genotypes-1.vcf"));
    }

    @BeforeClass
    public static void initializeTestDirectory() throws IOException {
        System.out.println("Creating base test directory: " + BASE_TEST_DIR);
        if (LOG.isDebugEnabled()) {
            LOG.debug("Creating base test directory: " + BASE_TEST_DIR);
        }
        FileUtils.forceMkdir(new File(BASE_TEST_DIR));
    }

    @Before
    public void setUp() throws IOException {

        if (LOG.isDebugEnabled()) {
            LOG.debug("Using test directory: " + BASE_TEST_DIR);
        }
    }


}
