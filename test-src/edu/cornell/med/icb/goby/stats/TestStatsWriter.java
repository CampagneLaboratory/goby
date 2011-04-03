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

        int fieldC=writer.defineField("INFO", "C", 1, ColumnType.String, "something about C");
        int fieldD=writer.defineField("INFO", "D", 1, ColumnType.String, "something about D");

        writer.writeHeader();
        writer.setInfo(fieldC,"1");
        writer.setInfo(fieldD,"dvalue");
        writer.writeRecord();
        writer.close();

        assertEquals(new File("test-data/stats-writer/1.vcf"), new File("test-results/stats-writer/1.vcf"));
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
