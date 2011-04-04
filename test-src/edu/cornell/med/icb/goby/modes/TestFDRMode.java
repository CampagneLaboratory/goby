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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.io.FileUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.io.File;

import com.martiansoftware.jsap.JSAPException;
import edu.cornell.med.icb.goby.util.TestFiles;

/**
 * @author Fabien Campagne
 *         Date: Apr 4, 2011
 *         Time: 11:15:54 AM
 */
public class TestFDRMode extends TestFiles {
    private static final Log LOG = LogFactory.getLog(TestFDRMode.class);
    private static final String BASE_TEST_DIR = "test-results/fdr-mode";

    @BeforeClass
    public static void initializeTestDirectory() throws IOException {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Creating base test directory: " + BASE_TEST_DIR);
        }

        FileUtils.forceMkdir(new File(BASE_TEST_DIR));
    }

    @Test
    public void mergeVCF1() throws IOException, JSAPException {

        FalseDiscoveryRateMode mode = new FalseDiscoveryRateMode();
        String[] args = {
                "--mode", "fdr",
                "--vcf",
                "test-data/fdr-mode/file1.vcf",
                "test-data/fdr-mode/file2.vcf",
                "test-data/fdr-mode/file3.vcf",
                "--output", "test-results/fdr-mode/combined-file.vcf",
        };
        mode.configure(args);
        mode.execute();
        assertEquals(new File("test-data/fdr-mode/expected-combined-1-2-3.vcf"), new File("test-results/fdr-mode/combined-file.vcf"));
    }

}
