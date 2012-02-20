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

import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import org.junit.Test;

import java.io.StringWriter;
import java.util.ArrayList;

import static org.junit.Assert.assertEquals;

/**
 * @author Nyasha Chambwe
 *         Date: 1/25/12
 *         Time: 12:06 PM
 */
public class TestAnnotationAveragingWriter {


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
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 1 result: ", "Chromosome\tStart\tEnd\tFeature\t" +
                "#C[sample1][CpG]\t#Cm[sample1][CpG]\t#Sites[sample1][CpG]\tMR[sample1][CpG]\t" +
                "#C[sample1][CpA]\t#Cm[sample1][CpA]\t#Sites[sample1][CpA]\tMR[sample1][CpA]\t" +
                "#C[sample1][CpC]\t#Cm[sample1][CpC]\t#Sites[sample1][CpC]\tMR[sample1][CpC]\t" +
                "#C[sample1][CpT]\t#Cm[sample1][CpT]\t#Sites[sample1][CpT]\tMR[sample1][CpT]\t" +
                "#C[sample1][CpN]\t#Cm[sample1][CpN]\t#Sites[sample1][CpN]\tMR[sample1][CpN]\n" +
                "Case1\t4\t8\tannotation0\t8\t17\t2\t68.0000\t" +
                "0\t0\t0\t\t" +
                "0\t0\t0\t\t" +
                "0\t0\t0\t\t" +
                "0\t0\t0\t\n", stringWriter.getBuffer().toString());
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
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 2 result: ",
                "Chromosome\tStart\tEnd\tFeature\t" +
                        "#C[sample1][CpG]\t#Cm[sample1][CpG]\t#Sites[sample1][CpG]\tMR[sample1][CpG]\t" +
                        "#C[sample1][CpA]\t#Cm[sample1][CpA]\t#Sites[sample1][CpA]\tMR[sample1][CpA]\t" +
                        "#C[sample1][CpC]\t#Cm[sample1][CpC]\t#Sites[sample1][CpC]\tMR[sample1][CpC]\t" +
                        "#C[sample1][CpT]\t#Cm[sample1][CpT]\t#Sites[sample1][CpT]\tMR[sample1][CpT]\t" +
                        "#C[sample1][CpN]\t#Cm[sample1][CpN]\t#Sites[sample1][CpN]\tMR[sample1][CpN]\n" +
                        "Case2\t5\t9\tannotation1\t" +
                        "8\t16\t2\t66.6667\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\n" +
                        "Case2\t13\t17\tannotation2\t" +
                        "17\t6\t2\t26.0870\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\n", stringWriter.getBuffer().toString());
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
        assertEquals("Test Case 3 result: ", "Chromosome\tStart\tEnd\tFeature\t" +
                "#C[sample1][CpG]\t#Cm[sample1][CpG]\t#Sites[sample1][CpG]\tMR[sample1][CpG]\t" +
                "#C[sample1][CpA]\t#Cm[sample1][CpA]\t#Sites[sample1][CpA]\tMR[sample1][CpA]\t" +
                "#C[sample1][CpC]\t#Cm[sample1][CpC]\t#Sites[sample1][CpC]\tMR[sample1][CpC]\t" +
                "#C[sample1][CpT]\t#Cm[sample1][CpT]\t#Sites[sample1][CpT]\tMR[sample1][CpT]\t" +
                "#C[sample1][CpN]\t#Cm[sample1][CpN]\t#Sites[sample1][CpN]\tMR[sample1][CpN]\n" +
                "Case3\t4\t21\tannotation4\t" +
                "17\t19\t3\t52.7778\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\n" +
                "Case3\t11\t26\tannotation6\t" +
                "43\t24\t5\t35.8209\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\n" +
                "Case3\t4\t30\tannotation3\t" +
                "57\t43\t8\t43.0000\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\n" +
                "Case3\t11\t30\tannotation5\t" +
                "23\t9\t3\t28.1250\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\n"
                , stringWriter.getBuffer().toString());
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
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 4 result: ",
                "Chromosome\tStart\tEnd\tFeature\t" +
                        "#C[sample1][CpG]\t#Cm[sample1][CpG]\t#Sites[sample1][CpG]\tMR[sample1][CpG]\t" +
                        "#C[sample2][CpG]\t#Cm[sample2][CpG]\t#Sites[sample2][CpG]\tMR[sample2][CpG]\t" +
                        "#C[sample3][CpG]\t#Cm[sample3][CpG]\t#Sites[sample3][CpG]\tMR[sample3][CpG]\t" +
                        "#C[sample1][CpA]\t#Cm[sample1][CpA]\t#Sites[sample1][CpA]\tMR[sample1][CpA]\t" +
                        "#C[sample2][CpA]\t#Cm[sample2][CpA]\t#Sites[sample2][CpA]\tMR[sample2][CpA]\t" +
                        "#C[sample3][CpA]\t#Cm[sample3][CpA]\t#Sites[sample3][CpA]\tMR[sample3][CpA]\t" +
                        "#C[sample1][CpC]\t#Cm[sample1][CpC]\t#Sites[sample1][CpC]\tMR[sample1][CpC]\t" +
                        "#C[sample2][CpC]\t#Cm[sample2][CpC]\t#Sites[sample2][CpC]\tMR[sample2][CpC]\t" +
                        "#C[sample3][CpC]\t#Cm[sample3][CpC]\t#Sites[sample3][CpC]\tMR[sample3][CpC]\t" +
                        "#C[sample1][CpT]\t#Cm[sample1][CpT]\t#Sites[sample1][CpT]\tMR[sample1][CpT]\t" +
                        "#C[sample2][CpT]\t#Cm[sample2][CpT]\t#Sites[sample2][CpT]\tMR[sample2][CpT]\t" +
                        "#C[sample3][CpT]\t#Cm[sample3][CpT]\t#Sites[sample3][CpT]\tMR[sample3][CpT]\t" +
                        "#C[sample1][CpN]\t#Cm[sample1][CpN]\t#Sites[sample1][CpN]\tMR[sample1][CpN]\t" +
                        "#C[sample2][CpN]\t#Cm[sample2][CpN]\t#Sites[sample2][CpN]\tMR[sample2][CpN]\t" +
                        "#C[sample3][CpN]\t#Cm[sample3][CpN]\t#Sites[sample3][CpN]\tMR[sample3][CpN]\n" +
                        "Case2\t5\t9\tannotation1\t" +
                        "8\t16\t2\t66.6667\t10\t16\t2\t61.5385\t11\t5\t2\t31.2500\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\n" +
                        "Case2\t13\t17\tannotation2\t" +
                        "17\t6\t2\t26.0870\t5\t12\t2\t70.5882\t17\t10\t2\t37.0370\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\n"
                , stringWriter.getBuffer().toString());
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
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        int[] a = {0, 0, 0};
        testWriter.setSampleIndexToGroupIndex(a);
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 5 result: ",
                "Chromosome\tStart\tEnd\tFeature\t" +
                        "#C[sample1][CpG]\t#Cm[sample1][CpG]\t#Sites[sample1][CpG]\tMR[sample1][CpG]\t" +
                        "#C[sample2][CpG]\t#Cm[sample2][CpG]\t#Sites[sample2][CpG]\tMR[sample2][CpG]\t" +
                        "#C[sample3][CpG]\t#Cm[sample3][CpG]\t#Sites[sample3][CpG]\tMR[sample3][CpG]\t" +
                        "#C[sample1][CpA]\t#Cm[sample1][CpA]\t#Sites[sample1][CpA]\tMR[sample1][CpA]\t" +
                        "#C[sample2][CpA]\t#Cm[sample2][CpA]\t#Sites[sample2][CpA]\tMR[sample2][CpA]\t" +
                        "#C[sample3][CpA]\t#Cm[sample3][CpA]\t#Sites[sample3][CpA]\tMR[sample3][CpA]\t" +
                        "#C[sample1][CpC]\t#Cm[sample1][CpC]\t#Sites[sample1][CpC]\tMR[sample1][CpC]\t" +
                        "#C[sample2][CpC]\t#Cm[sample2][CpC]\t#Sites[sample2][CpC]\tMR[sample2][CpC]\t" +
                        "#C[sample3][CpC]\t#Cm[sample3][CpC]\t#Sites[sample3][CpC]\tMR[sample3][CpC]\t" +
                        "#C[sample1][CpT]\t#Cm[sample1][CpT]\t#Sites[sample1][CpT]\tMR[sample1][CpT]\t" +
                        "#C[sample2][CpT]\t#Cm[sample2][CpT]\t#Sites[sample2][CpT]\tMR[sample2][CpT]\t" +
                        "#C[sample3][CpT]\t#Cm[sample3][CpT]\t#Sites[sample3][CpT]\tMR[sample3][CpT]\t" +
                        "#C[sample1][CpN]\t#Cm[sample1][CpN]\t#Sites[sample1][CpN]\tMR[sample1][CpN]\t" +
                        "#C[sample2][CpN]\t#Cm[sample2][CpN]\t#Sites[sample2][CpN]\tMR[sample2][CpN]\t" +
                        "#C[sample3][CpN]\t#Cm[sample3][CpN]\t#Sites[sample3][CpN]\tMR[sample3][CpN]\t" +
                        "#C[group1][CpG]\t#Cm[group1][CpG]\t#Sites[group1][CpG]\tMR[group1][CpG]\t" +
                        "#C[group1][CpA]\t#Cm[group1][CpA]\t#Sites[group1][CpA]\tMR[group1][CpA]\t" +
                        "#C[group1][CpC]\t#Cm[group1][CpC]\t#Sites[group1][CpC]\tMR[group1][CpC]\t" +
                        "#C[group1][CpT]\t#Cm[group1][CpT]\t#Sites[group1][CpT]\tMR[group1][CpT]\t" +
                        "#C[group1][CpN]\t#Cm[group1][CpN]\t#Sites[group1][CpN]\tMR[group1][CpN]\n" +
                        "Case2\t5\t9\tannotation1\t" +
                        "8\t16\t2\t66.6667\t" +
                        "10\t16\t2\t61.5385\t" +
                        "11\t5\t2\t31.2500\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "29\t37\t6\t56.0606\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\n" +
                        "Case2\t13\t17\tannotation2\t" +
                        "17\t6\t2\t26.0870\t" +
                        "5\t12\t2\t70.5882\t" +
                        "17\t10\t2\t37.0370\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "39\t28\t6\t41.7910\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\n"
                , stringWriter.getBuffer().toString());

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
                "Chromosome\tStart\tEnd\tFeature\t" +
                        "#C[sample1][CpG]\t#Cm[sample1][CpG]\t#Sites[sample1][CpG]\tMR[sample1][CpG]\t" +
                        "#C[sample2][CpG]\t#Cm[sample2][CpG]\t#Sites[sample2][CpG]\tMR[sample2][CpG]\t" +
                        "#C[sample3][CpG]\t#Cm[sample3][CpG]\t#Sites[sample3][CpG]\tMR[sample3][CpG]\t" +
                        "#C[sample4][CpG]\t#Cm[sample4][CpG]\t#Sites[sample4][CpG]\tMR[sample4][CpG]\t" +
                        "#C[sample5][CpG]\t#Cm[sample5][CpG]\t#Sites[sample5][CpG]\tMR[sample5][CpG]\t" +
                        "#C[sample6][CpG]\t#Cm[sample6][CpG]\t#Sites[sample6][CpG]\tMR[sample6][CpG]\t" +
                        "#C[sample1][CpA]\t#Cm[sample1][CpA]\t#Sites[sample1][CpA]\tMR[sample1][CpA]\t" +
                        "#C[sample2][CpA]\t#Cm[sample2][CpA]\t#Sites[sample2][CpA]\tMR[sample2][CpA]\t" +
                        "#C[sample3][CpA]\t#Cm[sample3][CpA]\t#Sites[sample3][CpA]\tMR[sample3][CpA]\t" +
                        "#C[sample4][CpA]\t#Cm[sample4][CpA]\t#Sites[sample4][CpA]\tMR[sample4][CpA]\t" +
                        "#C[sample5][CpA]\t#Cm[sample5][CpA]\t#Sites[sample5][CpA]\tMR[sample5][CpA]\t" +
                        "#C[sample6][CpA]\t#Cm[sample6][CpA]\t#Sites[sample6][CpA]\tMR[sample6][CpA]\t" +
                        "#C[sample1][CpC]\t#Cm[sample1][CpC]\t#Sites[sample1][CpC]\tMR[sample1][CpC]\t" +
                        "#C[sample2][CpC]\t#Cm[sample2][CpC]\t#Sites[sample2][CpC]\tMR[sample2][CpC]\t" +
                        "#C[sample3][CpC]\t#Cm[sample3][CpC]\t#Sites[sample3][CpC]\tMR[sample3][CpC]\t" +
                        "#C[sample4][CpC]\t#Cm[sample4][CpC]\t#Sites[sample4][CpC]\tMR[sample4][CpC]\t" +
                        "#C[sample5][CpC]\t#Cm[sample5][CpC]\t#Sites[sample5][CpC]\tMR[sample5][CpC]\t" +
                        "#C[sample6][CpC]\t#Cm[sample6][CpC]\t#Sites[sample6][CpC]\tMR[sample6][CpC]\t" +
                        "#C[sample1][CpT]\t#Cm[sample1][CpT]\t#Sites[sample1][CpT]\tMR[sample1][CpT]\t" +
                        "#C[sample2][CpT]\t#Cm[sample2][CpT]\t#Sites[sample2][CpT]\tMR[sample2][CpT]\t" +
                        "#C[sample3][CpT]\t#Cm[sample3][CpT]\t#Sites[sample3][CpT]\tMR[sample3][CpT]\t" +
                        "#C[sample4][CpT]\t#Cm[sample4][CpT]\t#Sites[sample4][CpT]\tMR[sample4][CpT]\t" +
                        "#C[sample5][CpT]\t#Cm[sample5][CpT]\t#Sites[sample5][CpT]\tMR[sample5][CpT]\t" +
                        "#C[sample6][CpT]\t#Cm[sample6][CpT]\t#Sites[sample6][CpT]\tMR[sample6][CpT]\t" +
                        "#C[sample1][CpN]\t#Cm[sample1][CpN]\t#Sites[sample1][CpN]\tMR[sample1][CpN]\t" +
                        "#C[sample2][CpN]\t#Cm[sample2][CpN]\t#Sites[sample2][CpN]\tMR[sample2][CpN]\t" +
                        "#C[sample3][CpN]\t#Cm[sample3][CpN]\t#Sites[sample3][CpN]\tMR[sample3][CpN]\t" +
                        "#C[sample4][CpN]\t#Cm[sample4][CpN]\t#Sites[sample4][CpN]\tMR[sample4][CpN]\t" +
                        "#C[sample5][CpN]\t#Cm[sample5][CpN]\t#Sites[sample5][CpN]\tMR[sample5][CpN]\t" +
                        "#C[sample6][CpN]\t#Cm[sample6][CpN]\t#Sites[sample6][CpN]\tMR[sample6][CpN]\t" +
                        "#C[group1][CpG]\t#Cm[group1][CpG]\t#Sites[group1][CpG]\tMR[group1][CpG]\t" +
                        "#C[group2][CpG]\t#Cm[group2][CpG]\t#Sites[group2][CpG]\tMR[group2][CpG]\t" +
                        "#C[group1][CpA]\t#Cm[group1][CpA]\t#Sites[group1][CpA]\tMR[group1][CpA]\t" +
                        "#C[group2][CpA]\t#Cm[group2][CpA]\t#Sites[group2][CpA]\tMR[group2][CpA]\t" +
                        "#C[group1][CpC]\t#Cm[group1][CpC]\t#Sites[group1][CpC]\tMR[group1][CpC]\t" +
                        "#C[group2][CpC]\t#Cm[group2][CpC]\t#Sites[group2][CpC]\tMR[group2][CpC]\t" +
                        "#C[group1][CpT]\t#Cm[group1][CpT]\t#Sites[group1][CpT]\tMR[group1][CpT]\t" +
                        "#C[group2][CpT]\t#Cm[group2][CpT]\t#Sites[group2][CpT]\tMR[group2][CpT]\t" +
                        "#C[group1][CpN]\t#Cm[group1][CpN]\t#Sites[group1][CpN]\tMR[group1][CpN]\t" +
                        "#C[group2][CpN]\t#Cm[group2][CpN]\t#Sites[group2][CpN]\tMR[group2][CpN]\n" +
                        "Case3\t4\t21\tannotation4\t" +
                        "13\t19\t3\t59.3750\t" +
                        "17\t16\t3\t48.4848\t" +
                        "15\t23\t3\t60.5263\t" +
                        "17\t12\t3\t41.3793\t" +
                        "10\t15\t3\t60.0000\t" +
                        "11\t18\t3\t62.0690\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "45\t58\t9\t56.3107\t" +
                        "38\t45\t9\t54.2169\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\n" +
                        "Case3\t11\t26\tannotation6\t" +
                        "27\t16\t5\t37.2093\t" +
                        "38\t19\t5\t33.3333\t" +
                        "41\t26\t5\t38.8060\t" +
                        "43\t27\t5\t38.5714\t" +
                        "36\t32\t5\t47.0588\t" +
                        "36\t24\t5\t40.0000\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "106\t61\t15\t36.5269\t" +
                        "115\t83\t15\t41.9192\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\n" +
                        "Case3\t4\t30\tannotation3\t" +
                        "42\t35\t8\t45.4545\t" +
                        "55\t34\t8\t38.2022\t" +
                        "49\t44\t8\t47.3118\t" +
                        "51\t40\t8\t43.9560\t" +
                        "44\t42\t8\t48.8372\t" +
                        "39\t42\t8\t51.8519\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "146\t113\t24\t43.6293\t" +
                        "134\t124\t24\t48.0620\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\n" +
                        "Case3\t11\t30\tannotation5\t" +
                        "19\t9\t3\t32.1429\t" +
                        "26\t11\t3\t29.7297\t" +
                        "17\t6\t3\t26.0870\t" +
                        "17\t13\t3\t43.3333\t" +
                        "17\t12\t3\t41.3793\t" +
                        "13\t9\t3\t40.9091\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "62\t26\t9\t29.5455\t47\t34\t9\t41.9753\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\n", stringWriter.getBuffer().toString());

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
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        int[] a = {0, 0, 0};
        testWriter.setSampleIndexToGroupIndex(a);
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 7 result: ",
                "Chromosome\tStart\tEnd\tFeature\t" +
                        "#C[sample1][CpG]\t#Cm[sample1][CpG]\t#Sites[sample1][CpG]\tMR[sample1][CpG]\t" +
                        "#C[sample2][CpG]\t#Cm[sample2][CpG]\t#Sites[sample2][CpG]\tMR[sample2][CpG]\t" +
                        "#C[sample3][CpG]\t#Cm[sample3][CpG]\t#Sites[sample3][CpG]\tMR[sample3][CpG]\t" +
                        "#C[sample1][CpA]\t#Cm[sample1][CpA]\t#Sites[sample1][CpA]\tMR[sample1][CpA]\t" +
                        "#C[sample2][CpA]\t#Cm[sample2][CpA]\t#Sites[sample2][CpA]\tMR[sample2][CpA]\t" +
                        "#C[sample3][CpA]\t#Cm[sample3][CpA]\t#Sites[sample3][CpA]\tMR[sample3][CpA]\t" +
                        "#C[sample1][CpC]\t#Cm[sample1][CpC]\t#Sites[sample1][CpC]\tMR[sample1][CpC]\t" +
                        "#C[sample2][CpC]\t#Cm[sample2][CpC]\t#Sites[sample2][CpC]\tMR[sample2][CpC]\t" +
                        "#C[sample3][CpC]\t#Cm[sample3][CpC]\t#Sites[sample3][CpC]\tMR[sample3][CpC]\t" +
                        "#C[sample1][CpT]\t#Cm[sample1][CpT]\t#Sites[sample1][CpT]\tMR[sample1][CpT]\t" +
                        "#C[sample2][CpT]\t#Cm[sample2][CpT]\t#Sites[sample2][CpT]\tMR[sample2][CpT]\t" +
                        "#C[sample3][CpT]\t#Cm[sample3][CpT]\t#Sites[sample3][CpT]\tMR[sample3][CpT]\t" +
                        "#C[sample1][CpN]\t#Cm[sample1][CpN]\t#Sites[sample1][CpN]\tMR[sample1][CpN]\t" +
                        "#C[sample2][CpN]\t#Cm[sample2][CpN]\t#Sites[sample2][CpN]\tMR[sample2][CpN]\t" +
                        "#C[sample3][CpN]\t#Cm[sample3][CpN]\t#Sites[sample3][CpN]\tMR[sample3][CpN]\t" +
                        "#C[group1][CpG]\t#Cm[group1][CpG]\t#Sites[group1][CpG]\tMR[group1][CpG]\t" +
                        "#C[group1][CpA]\t#Cm[group1][CpA]\t#Sites[group1][CpA]\tMR[group1][CpA]\t" +
                        "#C[group1][CpC]\t#Cm[group1][CpC]\t#Sites[group1][CpC]\tMR[group1][CpC]\t" +
                        "#C[group1][CpT]\t#Cm[group1][CpT]\t#Sites[group1][CpT]\tMR[group1][CpT]\t" +
                        "#C[group1][CpN]\t#Cm[group1][CpN]\t#Sites[group1][CpN]\tMR[group1][CpN]\n" +
                        "Case4\t5\t9\tannotation7\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "8\t16\t2\t66.6667\t" +
                        "10\t16\t2\t61.5385\t" +
                        "11\t5\t2\t31.2500\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "29\t37\t6\t56.0606\t" +
                        "0\t0\t0\t\n" +
                        "Case4\t13\t17\tannotation8\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "17\t6\t2\t26.0870\t" +
                        "5\t12\t2\t70.5882\t" +
                        "17\t10\t2\t37.0370\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +

                        "39\t28\t6\t41.7910\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\n"
                , stringWriter.getBuffer().toString());

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
                "Chromosome\tStart\tEnd\tFeature\t" +
                        "#C[sample1][CpG]\t#Cm[sample1][CpG]\t#Sites[sample1][CpG]\tMR[sample1][CpG]\t" +
                        "#C[sample2][CpG]\t#Cm[sample2][CpG]\t#Sites[sample2][CpG]\tMR[sample2][CpG]\t" +
                        "#C[sample3][CpG]\t#Cm[sample3][CpG]\t#Sites[sample3][CpG]\tMR[sample3][CpG]\t" +
                        "#C[sample4][CpG]\t#Cm[sample4][CpG]\t#Sites[sample4][CpG]\tMR[sample4][CpG]\t" +
                        "#C[sample1][CpA]\t#Cm[sample1][CpA]\t#Sites[sample1][CpA]\tMR[sample1][CpA]\t" +
                        "#C[sample2][CpA]\t#Cm[sample2][CpA]\t#Sites[sample2][CpA]\tMR[sample2][CpA]\t" +
                        "#C[sample3][CpA]\t#Cm[sample3][CpA]\t#Sites[sample3][CpA]\tMR[sample3][CpA]\t" +
                        "#C[sample4][CpA]\t#Cm[sample4][CpA]\t#Sites[sample4][CpA]\tMR[sample4][CpA]\t" +
                        "#C[sample1][CpC]\t#Cm[sample1][CpC]\t#Sites[sample1][CpC]\tMR[sample1][CpC]\t" +
                        "#C[sample2][CpC]\t#Cm[sample2][CpC]\t#Sites[sample2][CpC]\tMR[sample2][CpC]\t" +
                        "#C[sample3][CpC]\t#Cm[sample3][CpC]\t#Sites[sample3][CpC]\tMR[sample3][CpC]\t" +
                        "#C[sample4][CpC]\t#Cm[sample4][CpC]\t#Sites[sample4][CpC]\tMR[sample4][CpC]\t" +
                        "#C[sample1][CpT]\t#Cm[sample1][CpT]\t#Sites[sample1][CpT]\tMR[sample1][CpT]\t" +
                        "#C[sample2][CpT]\t#Cm[sample2][CpT]\t#Sites[sample2][CpT]\tMR[sample2][CpT]\t" +
                        "#C[sample3][CpT]\t#Cm[sample3][CpT]\t#Sites[sample3][CpT]\tMR[sample3][CpT]\t" +
                        "#C[sample4][CpT]\t#Cm[sample4][CpT]\t#Sites[sample4][CpT]\tMR[sample4][CpT]\t" +
                        "#C[sample1][CpN]\t#Cm[sample1][CpN]\t#Sites[sample1][CpN]\tMR[sample1][CpN]\t" +
                        "#C[sample2][CpN]\t#Cm[sample2][CpN]\t#Sites[sample2][CpN]\tMR[sample2][CpN]\t" +
                        "#C[sample3][CpN]\t#Cm[sample3][CpN]\t#Sites[sample3][CpN]\tMR[sample3][CpN]\t" +
                        "#C[sample4][CpN]\t#Cm[sample4][CpN]\t#Sites[sample4][CpN]\tMR[sample4][CpN]\t" +
                        "#C[group1][CpG]\t#Cm[group1][CpG]\t#Sites[group1][CpG]\tMR[group1][CpG]\t" +
                        "#C[group2][CpG]\t#Cm[group2][CpG]\t#Sites[group2][CpG]\tMR[group2][CpG]\t" +
                        "#C[group1][CpA]\t#Cm[group1][CpA]\t#Sites[group1][CpA]\tMR[group1][CpA]\t" +
                        "#C[group2][CpA]\t#Cm[group2][CpA]\t#Sites[group2][CpA]\tMR[group2][CpA]\t" +
                        "#C[group1][CpC]\t#Cm[group1][CpC]\t#Sites[group1][CpC]\tMR[group1][CpC]\t" +
                        "#C[group2][CpC]\t#Cm[group2][CpC]\t#Sites[group2][CpC]\tMR[group2][CpC]\t" +
                        "#C[group1][CpT]\t#Cm[group1][CpT]\t#Sites[group1][CpT]\tMR[group1][CpT]\t" +
                        "#C[group2][CpT]\t#Cm[group2][CpT]\t#Sites[group2][CpT]\tMR[group2][CpT]\t" +
                        "#C[group1][CpN]\t#Cm[group1][CpN]\t#Sites[group1][CpN]\tMR[group1][CpN]\t" +
                        "#C[group2][CpN]\t#Cm[group2][CpN]\t#Sites[group2][CpN]\tMR[group2][CpN]\t" +
                        "fisherP[group1/group2][CpG]\tfisherP[group1/group2][CpA]\t" +
                        "fisherP[group1/group2][CpC]\tfisherP[group1/group2][CpT]\t" +
                        "fisherP[group1/group2][CpN]\tdeltaMR[group1/group2][CpG]\t" +
                        "deltaMR[group1/group2][CpA]\tdeltaMR[group1/group2][CpC]\t" +
                        "deltaMR[group1/group2][CpT]\tdeltaMR[group1/group2][CpN]\n" +
                        "Case4\t5\t9\tannotation7\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "8\t16\t2\t66.6667\t" +
                        "10\t16\t2\t61.5385\t" +
                        "11\t5\t2\t31.2500\t" +
                        "11\t5\t2\t31.2500\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "18\t32\t4\t64.0000\t" +
                        "22\t10\t4\t31.2500\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "1.00000\t1.00000\t1.00000\t0.00622665\t1.00000\t\t\t\t32.7500\t\n" +
                        "Case4\t13\t17\tannotation8\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "17\t6\t2\t26.0870\t" +
                        "5\t12\t2\t70.5882\t" +
                        "17\t10\t2\t37.0370\t" +
                        "16\t9\t2\t36.0000\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "22\t18\t4\t45.0000\t" +
                        "33\t19\t4\t36.5385\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "0\t0\t0\t\t" +
                        "1.00000\t0.520511\t1.00000\t1.00000\t1.00000\t\t8.46154\t\t\t\n"
                , stringWriter.getBuffer().toString());
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
                        "#C[sample1][CpG]\t#Cm[sample1][CpG]\t#Sites[sample1][CpG]\tMR[sample1][CpG]\t" +
                        "#C[sample2][CpG]\t#Cm[sample2][CpG]\t#Sites[sample2][CpG]\tMR[sample2][CpG]\t" +
                        "#C[sample3][CpG]\t#Cm[sample3][CpG]\t#Sites[sample3][CpG]\tMR[sample3][CpG]\t" +
                        "#C[sample4][CpG]\t#Cm[sample4][CpG]\t#Sites[sample4][CpG]\tMR[sample4][CpG]\t" +

                        "#C[sample1][CpA]\t#Cm[sample1][CpA]\t#Sites[sample1][CpA]\tMR[sample1][CpA]\t" +
                        "#C[sample2][CpA]\t#Cm[sample2][CpA]\t#Sites[sample2][CpA]\tMR[sample2][CpA]\t" +
                        "#C[sample3][CpA]\t#Cm[sample3][CpA]\t#Sites[sample3][CpA]\tMR[sample3][CpA]\t" +
                        "#C[sample4][CpA]\t#Cm[sample4][CpA]\t#Sites[sample4][CpA]\tMR[sample4][CpA]\t" +

                        "#C[sample1][CpC]\t#Cm[sample1][CpC]\t#Sites[sample1][CpC]\tMR[sample1][CpC]\t" +
                        "#C[sample2][CpC]\t#Cm[sample2][CpC]\t#Sites[sample2][CpC]\tMR[sample2][CpC]\t" +
                        "#C[sample3][CpC]\t#Cm[sample3][CpC]\t#Sites[sample3][CpC]\tMR[sample3][CpC]\t" +
                        "#C[sample4][CpC]\t#Cm[sample4][CpC]\t#Sites[sample4][CpC]\tMR[sample4][CpC]\t" +

                        "#C[sample1][CpT]\t#Cm[sample1][CpT]\t#Sites[sample1][CpT]\tMR[sample1][CpT]\t" +
                        "#C[sample2][CpT]\t#Cm[sample2][CpT]\t#Sites[sample2][CpT]\tMR[sample2][CpT]\t" +
                        "#C[sample3][CpT]\t#Cm[sample3][CpT]\t#Sites[sample3][CpT]\tMR[sample3][CpT]\t" +
                        "#C[sample4][CpT]\t#Cm[sample4][CpT]\t#Sites[sample4][CpT]\tMR[sample4][CpT]\t" +

                        "#C[sample1][CpN]\t#Cm[sample1][CpN]\t#Sites[sample1][CpN]\tMR[sample1][CpN]\t" +
                        "#C[sample2][CpN]\t#Cm[sample2][CpN]\t#Sites[sample2][CpN]\tMR[sample2][CpN]\t" +
                        "#C[sample3][CpN]\t#Cm[sample3][CpN]\t#Sites[sample3][CpN]\tMR[sample3][CpN]\t" +
                        "#C[sample4][CpN]\t#Cm[sample4][CpN]\t#Sites[sample4][CpN]\tMR[sample4][CpN]\t" +

                        "#C[group1][CpG]\t#Cm[group1][CpG]\t#Sites[group1][CpG]\tMR[group1][CpG]\t" +
                        "#C[group2][CpG]\t#Cm[group2][CpG]\t#Sites[group2][CpG]\tMR[group2][CpG]\t" +
                        "#C[group3][CpG]\t#Cm[group3][CpG]\t#Sites[group3][CpG]\tMR[group3][CpG]\t" +
                        "#C[group4][CpG]\t#Cm[group4][CpG]\t#Sites[group4][CpG]\tMR[group4][CpG]\t" +

                        "#C[group1][CpA]\t#Cm[group1][CpA]\t#Sites[group1][CpA]\tMR[group1][CpA]\t" +
                        "#C[group2][CpA]\t#Cm[group2][CpA]\t#Sites[group2][CpA]\tMR[group2][CpA]\t" +
                        "#C[group3][CpA]\t#Cm[group3][CpA]\t#Sites[group3][CpA]\tMR[group3][CpA]\t" +
                        "#C[group4][CpA]\t#Cm[group4][CpA]\t#Sites[group4][CpA]\tMR[group4][CpA]\t" +

                        "#C[group1][CpC]\t#Cm[group1][CpC]\t#Sites[group1][CpC]\tMR[group1][CpC]\t" +
                        "#C[group2][CpC]\t#Cm[group2][CpC]\t#Sites[group2][CpC]\tMR[group2][CpC]\t" +
                        "#C[group3][CpC]\t#Cm[group3][CpC]\t#Sites[group3][CpC]\tMR[group3][CpC]\t" +
                        "#C[group4][CpC]\t#Cm[group4][CpC]\t#Sites[group4][CpC]\tMR[group4][CpC]\t" +

                        "#C[group1][CpT]\t#Cm[group1][CpT]\t#Sites[group1][CpT]\tMR[group1][CpT]\t" +
                        "#C[group2][CpT]\t#Cm[group2][CpT]\t#Sites[group2][CpT]\tMR[group2][CpT]\t" +
                        "#C[group3][CpT]\t#Cm[group3][CpT]\t#Sites[group3][CpT]\tMR[group3][CpT]\t" +
                        "#C[group4][CpT]\t#Cm[group4][CpT]\t#Sites[group4][CpT]\tMR[group4][CpT]\t" +

                        "#C[group1][CpN]\t#Cm[group1][CpN]\t#Sites[group1][CpN]\tMR[group1][CpN]\t" +
                        "#C[group2][CpN]\t#Cm[group2][CpN]\t#Sites[group2][CpN]\tMR[group2][CpN]\t" +
                        "#C[group3][CpN]\t#Cm[group3][CpN]\t#Sites[group3][CpN]\tMR[group3][CpN]\t" +
                        "#C[group4][CpN]\t#Cm[group4][CpN]\t#Sites[group4][CpN]\tMR[group4][CpN]\t" +

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

                        "Case4\t5\t9\tannotation7\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "8\t16\t2\t66.6667\t" +
                        "10\t16\t2\t61.5385\t" +
                        "11\t5\t2\t31.2500\t" +
                        "11\t5\t2\t31.2500\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "8\t16\t2\t66.6667\t" +
                        "10\t16\t2\t61.5385\t" +
                        "11\t5\t2\t31.2500\t" +
                        "11\t5\t2\t31.2500\t"+
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t0.773708\t0.0514794\t0.0514794\t1.00000\t1.00000\t1.00000\t\t\t\t\t\t\t\t\t\t5.12821\t35.4167\t35.4167\t\t\t\n" +

                        "Case4\t13\t17\tannotation8\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "17\t6\t2\t26.0870\t5\t12\t2\t70.5882\t17\t10\t2\t37.0370\t16\t9\t2\t36.0000\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "17\t6\t2\t26.0870\t5\t12\t2\t70.5882\t17\t10\t2\t37.0370\t16\t9\t2\t36.0000\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0\t\t" +
                        "1.00000\t1.00000\t1.00000\t0.00952120\t0.545556\t0.541875\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t1.00000\t\t\t\t44.5013\t10.9501\t9.91304\t\t\t\t\t\t\t\t\t\n"
                ,stringWriter.getBuffer().toString());

    }
}
