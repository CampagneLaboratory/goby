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

package edu.cornell.med.icb.goby.stats;

import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.EstimatedDistribution;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.FastSmallAndLog10BinningStrategy;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.PassThroughStatisticAdaptor;
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;

import static org.junit.Assert.assertEquals;

/**
 * @author Nyasha Chambwe
 *         Date: 1/25/12
 *         Time: 12:06 PM
 */
public class TestAnnotationAveragingWriter {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestAnnotationAveragingWriter.class);

    private static final String BASE_TEST_DIR = "test-results/stats";

    String[] sequences = {
            "GTTACGCGATGATTTAAGAAT",
            "TTAATCGCGTAATCGCGATAAT",
            "TAATCGCGTAACGCGATACGATTCGCGACGT",
            "TTAATCTCTTAATCACATAAT",
    };

    String[] sequenceNames = {"Case1", "Case2", "Case3", "Case4"};

    RandomAccessSequenceTestSupport genome = new RandomAccessSequenceTestSupport(sequences) {
        @Override
        public int getReferenceIndex(String referenceId) {

            if (referenceId.equals("Case1")) {
                return 0;
            }
            if (referenceId.equals("Case2")) {
                return 1;
            }
            if (referenceId.equals("Case3")) {
                return 2;
            }
            if (referenceId.equals("Case4")) {
                return 3;
            }

            return -1;
        }

        @Override
        public String getReferenceName(int index) {
            return sequenceNames[index];
        }
    };


    private MethylCountProvider testSupport;
    private static final String[] DEFAULT_TEST_CONTEXTS = {"CpG", "CpA", "CpC", "CpT", "CpN"};


    @Test
    public void testCase1() {
        String[] samples = new String[]{"sample1"};
        int[] positions = new int[]{5, 7};
        int[][] C = {{5, 3}};
        int[][] Cm = {{9, 8}};
        testSupport = new MethylCountProviderTestSupport(samples, positions, "Case1", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        testWriter.setContexts(DEFAULT_TEST_CONTEXTS);
        testWriter.setWriteNumSites(false);
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 1 result: ", "Chromosome\tStart\tEnd\tFeature\tMR[sample1][CpG]\tMR[sample1][CpA]\tMR[sample1][CpC]\tMR[sample1][CpT]\tMR[sample1][CpN]\n" +
                "Case1\t4\t8\tannotation0\t68.00\t\t\t\t\n", stringWriter.getBuffer().toString());
    }

    @Test
    public void testCase2() {
        String[] samples = new String[]{"sample1"};
        int[] positions = new int[]{6, 8, 14, 16};
        int[][] C = {{5, 3, 9, 8}};
        int[][] Cm = {{9, 7, 1, 5}};

        testSupport = new MethylCountProviderTestSupport(samples, positions, "Case2", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        testWriter.setWriteNumSites(false);
        testWriter.setContexts(DEFAULT_TEST_CONTEXTS);
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 2 result: ",
                "Chromosome\tStart\tEnd\tFeature\tMR[sample1][CpG]\tMR[sample1][CpA]\tMR[sample1][CpC]\tMR[sample1][CpT]\tMR[sample1][CpN]\n" +
                        "Case2\t5\t9\tannotation1\t66.67\t\t\t\t\n" +
                        "Case2\t13\t17\tannotation2\t26.09\t\t\t\t\n", stringWriter.getBuffer().toString());
    }

    @Test
    public void testCase3() {
        String[] samples = new String[]{"sample1"};
        int[] positions = new int[]{5, 7, 12, 14, 19, 24, 26, 29};
        int[][] C = {{5, 3, 9, 8, 9, 9, 8, 6}};
        int[][] Cm = {{9, 7, 1, 5, 3, 7, 8, 3}};
        testSupport = new MethylCountProviderTestSupport(samples, positions, "Case3", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setWriteNumSites(false);
        testWriter.setContexts(DEFAULT_TEST_CONTEXTS);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 3 result: ", "Chromosome\tStart\tEnd\tFeature\tMR[sample1][CpG]\tMR[sample1][CpA]\tMR[sample1][CpC]\tMR[sample1][CpT]\tMR[sample1][CpN]\n" +
                "Case3\t4\t21\tannotation4\t52.78\t\t\t\t\n" +
                "Case3\t11\t26\tannotation6\t35.82\t\t\t\t\n" +
                "Case3\t4\t30\tannotation3\t43.00\t\t\t\t\n" +
                "Case3\t11\t30\tannotation5\t28.13\t\t\t\t\n", stringWriter.getBuffer().toString());
    }

    @Test
    public void testCase4() {
        String[] samples = new String[]{"sample1", "sample2", "sample3"};
        int[] positions = new int[]{6, 8, 14, 16};
        int[][] C = {{5, 3, 9, 8}, {4, 6, 3, 2}, {8, 3, 8, 9}};
        int[][] Cm = {{9, 7, 1, 5}, {9, 7, 9, 3}, {2, 3, 2, 8}};
        testSupport = new MethylCountProviderTestSupport(samples, positions, "Case2", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setWriteNumSites(false);
        testWriter.setContexts(DEFAULT_TEST_CONTEXTS);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 4 result: ",
                "Chromosome\tStart\tEnd\tFeature\tMR[sample1][CpG]\tMR[sample2][CpG]\tMR[sample3][CpG]\tMR[sample1][CpA]\tMR[sample2][CpA]\tMR[sample3][CpA]\tMR[sample1][CpC]\tMR[sample2][CpC]\tMR[sample3][CpC]\tMR[sample1][CpT]\tMR[sample2][CpT]\tMR[sample3][CpT]\tMR[sample1][CpN]\tMR[sample2][CpN]\tMR[sample3][CpN]\n" +
                        "Case2\t5\t9\tannotation1\t66.67\t61.54\t31.25\t\t\t\t\t\t\t\t\t\t\t\t\n" +
                        "Case2\t13\t17\tannotation2\t26.09\t70.59\t37.04\n", stringWriter.getBuffer().toString());
    }

    @Test
    public void testCase5() {
        String[] groups = new String[]{"group1"};
        String[] samples = new String[]{"sample1", "sample2", "sample3"};
        int[] positions = new int[]{6, 8, 14, 16};
        int[][] C = {{5, 3, 9, 8}, {4, 6, 3, 2}, {8, 3, 8, 9}};
        int[][] Cm = {{9, 7, 1, 5}, {9, 7, 9, 3}, {2, 3, 2, 8}};
        testSupport = new MethylCountProviderTestSupport(groups, samples, positions, "Case2", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setWriteNumSites(false);
        testWriter.setContexts(DEFAULT_TEST_CONTEXTS);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        int[] a = {0, 0, 0};
        testWriter.setSampleIndexToGroupIndex(a);
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 5 result: ",
                "Chromosome\tStart\tEnd\tFeature\tMR[sample1][CpG]\tMR[sample2][CpG]\tMR[sample3][CpG]\tMR[sample1][CpA]\tMR[sample2][CpA]\tMR[sample3][CpA]\tMR[sample1][CpC]\tMR[sample2][CpC]\tMR[sample3][CpC]\tMR[sample1][CpT]\tMR[sample2][CpT]\tMR[sample3][CpT]\tMR[sample1][CpN]\tMR[sample2][CpN]\tMR[sample3][CpN]\tMR[group1][CpG]\tMR[group1][CpA]\tMR[group1][CpC]\tMR[group1][CpT]\tMR[group1][CpN]\n" +
                        "Case2\t5\t9\tannotation1\t66.67\t61.54\t31.25\t\t\t\t\t\t\t\t\t\t\t\t\t56.06\t\t\t\t\n" +
                        "Case2\t13\t17\tannotation2\t26.09\t70.59\t37.04\t\t\t\t\t\t\t\t\t\t\t\t\t41.79\t\t\t\tn", stringWriter.getBuffer().toString());

    }

    @Test
    public void testCase6() {
        String[] groups = new String[]{"group1", "group2"};
        String[] samples = new String[]{"sample1", "sample2", "sample3", "sample4", "sample5", "sample6"};
        int[] positions = new int[]{5, 7, 12, 14, 19, 24, 26, 29};
        int[][] C = {{5, 3, 4, 8, 5, 2, 8, 7}, {5, 3, 9, 8, 9, 8, 4, 9}, {5, 3, 9, 8, 7, 9, 8, 0}, {5, 3, 9, 8, 9, 9, 8, 0},
                {5, 3, 9, 8, 2, 9, 8, 0}, {0, 3, 5, 8, 8, 7, 8, 0}};
        int[][] Cm = {{9, 7, 1, 5, 3, 7, 0, 3}, {9, 3, 1, 7, 4, 7, 0, 3}, {9, 6, 1, 2, 8, 7, 8, 3}, {9, 1, 1, 9, 2, 7, 8, 3},
                {0, 7, 1, 8, 8, 7, 8, 3}, {8, 7, 1, 5, 3, 7, 8, 3}};
        testSupport = new MethylCountProviderTestSupport(groups, samples, positions, "Case3", C, Cm);

        final StringWriter stringWriter = new StringWriter();
        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setWriteNumSites(false);
        testWriter.setContexts(DEFAULT_TEST_CONTEXTS);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        int[] a = {0, 0, 0, 1, 1, 1};
        testWriter.setSampleIndexToGroupIndex(a);

        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 6 result: ",
                "Chromosome\tStart\tEnd\tFeature\tMR[sample1][CpG]\tMR[sample2][CpG]\tMR[sample3][CpG]\tMR[sample4][CpG]\tMR[sample5][CpG]\tMR[sample6][CpG]\tMR[sample1][CpA]\tMR[sample2][CpA]\tMR[sample3][CpA]\tMR[sample4][CpA]\tMR[sample5][CpA]\tMR[sample6][CpA]\tMR[sample1][CpC]\tMR[sample2][CpC]\tMR[sample3][CpC]\tMR[sample4][CpC]\tMR[sample5][CpC]\tMR[sample6][CpC]\tMR[sample1][CpT]\tMR[sample2][CpT]\tMR[sample3][CpT]\tMR[sample4][CpT]\tMR[sample5][CpT]\tMR[sample6][CpT]\tMR[sample1][CpN]\tMR[sample2][CpN]\tMR[sample3][CpN]\tMR[sample4][CpN]\tMR[sample5][CpN]\tMR[sample6][CpN]\tMR[group1][CpG]\tMR[group2][CpG]\tMR[group1][CpA]\tMR[group2][CpA]\tMR[group1][CpC]\tMR[group2][CpC]\tMR[group1][CpT]\tMR[group2][CpT]\tMR[group1][CpN]\tMR[group2][CpN]\n" +
                        "Case3\t4\t21\tannotation4\t59.38\t48.48\t60.53\t41.38\t60.00\t62.07\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t56.31\t54.22\t\t\t\t\t\t\t\t\n" +
                        "Case3\t11\t26\tannotation6\t37.21\t33.33\t38.81\t38.57\t47.06\t40.00\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t36.53\t41.92\t\t\t\t\t\t\t\t\n" +
                        "Case3\t4\t30\tannotation3\t45.45\t38.20\t47.31\t43.96\t48.84\t51.85\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t43.63\t48.06\t\t\t\t\t\t\t\t\n" +
                        "Case3\t11\t30\tannotation5\t32.14\t29.73\t26.09\t43.33\t41.38\t40.91\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t29.55\t41.98\t\t\t\t\t\t\t\tn", stringWriter.getBuffer().toString());

    }


    @Test
    public void testCase7() {
        String[] groups = new String[]{"group1"};
        String[] samples = new String[]{"sample1", "sample2", "sample3"};
        int[] positions = new int[]{6, 8, 14, 16};
        int[][] C = {{5, 3, 9, 8}, {4, 6, 3, 2}, {8, 3, 8, 9}};
        int[][] Cm = {{9, 7, 1, 5}, {9, 7, 9, 3}, {2, 3, 2, 8}};
        testSupport = new MethylCountProviderTestSupport(groups, samples, positions, "Case4", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setWriteNumSites(false);
        testWriter.setContexts(DEFAULT_TEST_CONTEXTS);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        int[] a = {0, 0, 0};
        testWriter.setSampleIndexToGroupIndex(a);
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 7 result: ",
                "Chromosome\tStart\tEnd\tFeature\tMR[sample1][CpG]\tMR[sample2][CpG]\tMR[sample3][CpG]\tMR[sample1][CpA]\tMR[sample2][CpA]\tMR[sample3][CpA]\tMR[sample1][CpC]\tMR[sample2][CpC]\tMR[sample3][CpC]\tMR[sample1][CpT]\tMR[sample2][CpT]\tMR[sample3][CpT]\tMR[sample1][CpN]\tMR[sample2][CpN]\tMR[sample3][CpN]\tMR[group1][CpG]\tMR[group1][CpA]\tMR[group1][CpC]\tMR[group1][CpT]\tMR[group1][CpN]\n" +
                        "Case4\t5\t9\tannotation7\t\t\t\t\t\t\t\t\t\t66.67\t61.54\t31.25\t\t\t\t\t\t\t56.06\t\n" +
                        "Case4\t13\t17\tannotation8\t\t\t\t26.09\t70.59\t37.04\t\t\t\t\t\t\t\t\t\t\t41.79\t\t\t\n", stringWriter.getBuffer().toString());

    }

    @Test
    public void testCase8() {
        // Test multiple samples, multiple groups, one comparison
        String[] groups = new String[]{"group1", "group2"};
        String[] samples = new String[]{"sample1", "sample2", "sample3", "sample4"};
        int[] positions = new int[]{6, 8, 14, 16};
        int[][] C = {{5, 3, 9, 8}, {4, 6, 3, 2}, {8, 3, 8, 9}, {7, 4, 9, 7}};
        int[][] Cm = {{9, 7, 1, 5}, {9, 7, 9, 3}, {2, 3, 2, 8}, {4, 1, 3, 6}};
        testSupport = new MethylCountProviderTestSupport(groups, samples, positions, "Case4", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setWriteNumSites(false);
        testWriter.setContexts(DEFAULT_TEST_CONTEXTS);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        int[] a = {0, 0, 1, 1};
        testWriter.setSampleIndexToGroupIndex(a);
        ArrayList<GroupComparison> groupComparisons = new ArrayList<GroupComparison>();
        GroupComparison case8 = new GroupComparison("group1", "group2", 0, 1, 0);
        groupComparisons.add(case8);
        testWriter.setGroupComparisons(groupComparisons);
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 8 result: ",
                "Chromosome\tStart\tEnd\tFeature\tMR[sample1][CpG]\tMR[sample2][CpG]\tMR[sample3][CpG]\tMR[sample4][CpG]\tMR[sample1][CpA]\tMR[sample2][CpA]\tMR[sample3][CpA]\tMR[sample4][CpA]\tMR[sample1][CpC]\tMR[sample2][CpC]\tMR[sample3][CpC]\tMR[sample4][CpC]\tMR[sample1][CpT]\tMR[sample2][CpT]\tMR[sample3][CpT]\tMR[sample4][CpT]\tMR[sample1][CpN]\tMR[sample2][CpN]\tMR[sample3][CpN]\tMR[sample4][CpN]\tMR[group1][CpG]\tMR[group2][CpG]\tMR[group1][CpA]\tMR[group2][CpA]\tMR[group1][CpC]\tMR[group2][CpC]\tMR[group1][CpT]\tMR[group2][CpT]\tMR[group1][CpN]\tMR[group2][CpN]\tfisherP[group1/group2][CpG]\tfisherP[group1/group2][CpA]\tfisherP[group1/group2][CpC]\tfisherP[group1/group2][CpT]\tfisherP[group1/group2][CpN]\tdeltaMR[group1/group2][CpG]\tdeltaMR[group1/group2][CpA]\tdeltaMR[group1/group2][CpC]\tdeltaMR[group1/group2][CpT]\tdeltaMR[group1/group2][CpN]\n" +
                        "Case4\t5\t9\tannotation7\t\t\t\t\t\t\t\t\t\t\t\t\t66.67\t61.54\t31.25\t31.25\t\t\t\t\t\t\t\t\t\t\t64.00\t31.25\t\t\t1.000\t1.000\t1.000\t0.006227\t1.000\t\t\t\t32.75\t\n" +
                        "Case4\t13\t17\tannotation8\t\t\t\t\t26.09\t70.59\t37.04\t36.00\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t45.00\t36.54\t\t\t\t\t\t\t1.000\t0.5205\t1.000\t1.000\t1.000\t\t8.462\t\t\t\n", stringWriter.getBuffer().toString());
    }

    @Test
    public void testCase9() {
        // test multiple samples, 4 groups, 3 comparisons
        // Test multiple samples, multiple groups, one comparison
        String[] groups = new String[]{"group1", "group2", "group3", "group4"};
        String[] samples = new String[]{"sample1", "sample2", "sample3", "sample4"};
        int[] positions = new int[]{6, 8, 14, 16};
        int[][] C = {{5, 3, 9, 8}, {4, 6, 3, 2}, {8, 3, 8, 9}, {7, 4, 9, 7}};
        int[][] Cm = {{9, 7, 1, 5}, {9, 7, 9, 3}, {2, 3, 2, 8}, {4, 1, 3, 6}};
        testSupport = new MethylCountProviderTestSupport(groups, samples, positions, "Case4", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setWriteNumSites(false);
        testWriter.setContexts(DEFAULT_TEST_CONTEXTS);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        int[] a = {0, 1, 2, 3};
        testWriter.setSampleIndexToGroupIndex(a);
        ArrayList<GroupComparison> groupComparisons = new ArrayList<GroupComparison>();
        GroupComparison comparisonToMake = new GroupComparison("group1", "group2", 0, 1, 0);
        groupComparisons.add(comparisonToMake);
        comparisonToMake = new GroupComparison("group1", "group3", 0, 2, 1);
        groupComparisons.add(comparisonToMake);
        comparisonToMake = new GroupComparison("group1", "group4", 0, 3, 2);
        groupComparisons.add(comparisonToMake);
        testWriter.setGroupComparisons(groupComparisons);
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 9 result: ",
                "Chromosome\tStart\tEnd\tFeature\t" +
                        "MR[sample1][CpG]\tMR[sample2][CpG]\tMR[sample3][CpG]\tMR[sample4][CpG]\t" +
                        "MR[sample1][CpA]\tMR[sample2][CpA]\tMR[sample3][CpA]\tMR[sample4][CpA]\t" +
                        "MR[sample1][CpC]\tMR[sample2][CpC]\tMR[sample3][CpC]\tMR[sample4][CpC]\t" +
                        "MR[sample1][CpT]\tMR[sample2][CpT]\tMR[sample3][CpT]\tMR[sample4][CpT]\t" +
                        "MR[sample1][CpN]\tMR[sample2][CpN]\tMR[sample3][CpN]\tMR[sample4][CpN]\t" +
                        "MR[group1][CpG]\tMR[group2][CpG]\tMR[group3][CpG]\tMR[group4][CpG]\t" +
                        "MR[group1][CpA]\tMR[group2][CpA]\tMR[group3][CpA]\tMR[group4][CpA]\t" +
                        "MR[group1][CpC]\tMR[group2][CpC]\tMR[group3][CpC]\tMR[group4][CpC]\t" +
                        "MR[group1][CpT]\tMR[group2][CpT]\tMR[group3][CpT]\tMR[group4][CpT]\t" +
                        "MR[group1][CpN]\tMR[group2][CpN]\tMR[group3][CpN]\tMR[group4][CpN]\t" +
                        "fisherP[group1/group2][CpG]\tfisherP[group1/group3][CpG]\tfisherP[group1/group4][CpG]\t" +
                        "fisherP[group1/group2][CpA]\tfisherP[group1/group3][CpA]\tfisherP[group1/group4][CpA]\t" +
                        "fisherP[group1/group2][CpC]\tfisherP[group1/group3][CpC]\tfisherP[group1/group4][CpC]\t" +
                        "fisherP[group1/group2][CpT]\tfisherP[group1/group3][CpT]\tfisherP[group1/group4][CpT]\t" +
                        "fisherP[group1/group2][CpN]\tfisherP[group1/group3][CpN]\tfisherP[group1/group4][CpN]\t" +
                        "deltaMR[group1/group2][CpG]\tdeltaMR[group1/group3][CpG]\tdeltaMR[group1/group4][CpG]\t" +
                        "deltaMR[group1/group2][CpA]\tdeltaMR[group1/group3][CpA]\tdeltaMR[group1/group4][CpA]\t" +
                        "deltaMR[group1/group2][CpC]\tdeltaMR[group1/group3][CpC]\tdeltaMR[group1/group4][CpC]\t" +
                        "deltaMR[group1/group2][CpT]\tdeltaMR[group1/group3][CpT]\tdeltaMR[group1/group4][CpT]\t" +
                        "deltaMR[group1/group2][CpN]\tdeltaMR[group1/group3][CpN]\tdeltaMR[group1/group4][CpN]\n" +
                        "Case4\t5\t9\tannotation7\t\t\t\t\t\t\t\t\t\t\t\t\t66.6667\t61.5385\t31.2500\t31.2500\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t66.6667\t61.5385\t31.2500\t31.2500\t\t\t\t\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t0.773708\t0.0514794\t0.0514794\t1.00000\t1.00000\t1.00000\t\t\t\t\t\t\t\t\t\t5.12821\t35.4167\t35.4167\t\t\t\n" +
                        "Case4\t13\t17\tannotation8\t\t\t\t\t26.0870\t70.5882\t37.0370\t36.0000\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t26.0870\t70.5882\t37.0370\t36.0000\t\t\t\t\t\t\t\t\t\t\t\t\t1.00000\t1.00000\t1.00000\t0.00952120\t0.545556\t0.541875\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t\t\t\t44.5013\t10.9501\t9.91304\t\t\t\t\t\t\t\t\t\n",
                stringWriter.getBuffer().toString());

    }

    @Test
    // test writing counts and sites
    public void testCase10() {
        String[] samples = new String[]{"sample1"};
        int[] positions = new int[]{5, 7};
        int[][] C = {{5, 3}};
        int[][] Cm = {{9, 8}};
        testSupport = new MethylCountProviderTestSupport(samples, positions, "Case1", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setWriteNumSites(true);
        testWriter.setWriteCounts(true);
        testWriter.setContexts(DEFAULT_TEST_CONTEXTS);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 1 result: ",
                "Chromosome\tStart\tEnd\tFeature\t" +
                        "#C[sample1][CpG]\t#C[sample1][CpA]\t#C[sample1][CpC]\t#C[sample1][CpT]\t#C[sample1][CpN]\t" +
                        "#Cm[sample1][CpG]\t#Cm[sample1][CpA]\t#Cm[sample1][CpC]\t#Cm[sample1][CpT]\t#Cm[sample1][CpN]\t" +
                        "MR[sample1][CpG]\tMR[sample1][CpA]\tMR[sample1][CpC]\tMR[sample1][CpT]\tMR[sample1][CpN]\t" +
                        "#Sites[sample1][CpG]\t#Sites[sample1][CpA]\t#Sites[sample1][CpC]\t#Sites[sample1][CpT]\t#Sites[sample1][CpN]\n" +
                        "Case1\t4\t8\tannotation0\t8\t0\t0\t0\t0\t17\t0\t0\t0\t0\t68.0000\t\t\t\t\t2\t0\t0\t0\t0\n", stringWriter.getBuffer().toString());
    }

    @Test
    // test writing counts
    public void testCaseEstimateP() throws IOException {

        String[] groups = new String[]{"group1", "group1", "group2", "group2"};
        String[] samples = new String[]{"sample1", "sample2", "sample3", "sample4"};
        int[] positions = new int[]{6, 8, 14, 16};
        int[][] C = {{5, 3, 9, 8}, {4, 6, 3, 2}, {8, 3, 8, 9}, {7, 4, 9, 7}};
        int[][] Cm = {{9, 7, 1, 5}, {9, 7, 9, 3}, {2, 3, 2, 8}, {4, 1, 3, 6}};
        testSupport = new MethylCountProviderTestSupport(groups, samples, positions, "Case4", C, Cm);
        final StringWriter stringWriter = new StringWriter();

        AnnotationAveragingWriter.doc().setValue("estimate-empirical-P", true);
        AnnotationAveragingWriter.doc().setValue("serialized-estimator-filename", makeDensityEstimator());
        AnnotationAveragingWriter.doc().setValue("combinator", "sum");

        AnnotationAveragingWriter testWriter = new AnnotationAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setWriteNumSites(false);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        testWriter.setAggregateAllContexts(true);
        testWriter.setSampleIndexToGroupIndex(new int[]{0, 0, 1, 1});
        ArrayList<GroupComparison> groupComparisons = new ArrayList<GroupComparison>();
        GroupComparison comparisonToMake = new GroupComparison("group1", "group2", 0, 1, 0);
        groupComparisons.add(comparisonToMake);
        testWriter.setGroupComparisons(groupComparisons);
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 1 result: ",
                "Chromosome\tStart\tEnd\tFeature\tMR[sample1]\tMR[sample2]\tMR[sample3]\tMR[sample4]\tMR[group1]\tMR[group1]\tMR[group2]\tMR[group2]\tfisherP[group1/group2]\tdeltaMR[group1/group2]\tempiricalP[group1/group2]\n" +
                        "Case4\t5\t9\tannotation7\t66.6667\t61.5385\t31.2500\t31.2500\t64.0000\t31.2500\t\t\t0.00622665\t32.7500\t0.000999500\n" +
                        "Case4\t13\t17\tannotation8\t26.0870\t70.5882\t37.0370\t36.0000\t45.0000\t36.5385\t\t\t0.520511\t8.46154\t0.000999500\n", stringWriter.getBuffer().toString());
    }

    private String makeDensityEstimator() throws IOException {
        int fixedValue = 40;
        EstimatedDistribution estimator = new EstimatedDistribution(1, new PassThroughStatisticAdaptor(100, fixedValue));
        estimator.setBinningStrategy(new FastSmallAndLog10BinningStrategy());
        int[] delta = {0000, 001, 002, 003, 004, 5};
        int[] frequencies = {1000, 900, 800, 700, 600, 1};
        for (int index = 0; index < delta.length; index++) {

            for (int i = 0; i < frequencies[index]; i++) {
                estimator.observe(0, fixedValue, delta[index]);
            }
        }
        String filename = FilenameUtils.concat(BASE_TEST_DIR, "density.bin");
        EstimatedDistribution.store(estimator, filename);
        return filename;
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
