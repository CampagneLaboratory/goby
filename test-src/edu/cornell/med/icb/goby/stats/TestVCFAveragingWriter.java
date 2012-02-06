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

import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import org.junit.Test;

import java.io.StringWriter;

import static org.junit.Assert.assertEquals;

/**
 * @author Nyasha Chambwe
 *         Date: 1/25/12
 *         Time: 12:06 PM
 */
public class TestVCFAveragingWriter {


    String[] sequences = {
            "Case1",
            "Case2",
            "Case3"
    };

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

            return -1;
        }

        @Override
        public String getReferenceName(int index) {
            return sequences[index];
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
        VCFAveragingWriter testWriter = new VCFAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 1 result: ", "Chromosome\tStart\tEnd\tFeature\tMR[sample1]\n" +
                "Case1\t4\t8\tannotation0\t68.0000\n", stringWriter.getBuffer().toString());
    }

    @Test
    public void testCase2() {
        String[] samples = new String[]{"sample1"};
        int[] positions = new int[]{6, 8, 14, 16};
        int[][] C = {{5, 3, 9, 8}};
        int[][] Cm = {{9, 7, 1, 5}};

        testSupport = new MethylCountProviderTestSupport(samples, positions, "Case2", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        VCFAveragingWriter testWriter = new VCFAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 2 result: ",
                "Chromosome\tStart\tEnd\tFeature\tMR[sample1]\n" +
                        "Case2\t5\t9\tannotation1\t66.6667\n" +
                        "Case2\t13\t17\tannotation2\t26.0870\n", stringWriter.getBuffer().toString());
    }

    @Test
    public void testCase3() {
        String[] samples = new String[]{"sample1"};
        int[] positions = new int[]{5, 7, 12, 14, 19, 24, 26, 29};
        int[][] C = {{5, 3, 9, 8, 9, 9, 8, 6}};
        int[][] Cm = {{9, 7, 1, 5, 3, 7, 8, 3}};
        testSupport = new MethylCountProviderTestSupport(samples, positions, "Case3", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        VCFAveragingWriter testWriter = new VCFAveragingWriter(stringWriter, genome, testSupport);
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

        assertEquals("Test Case 3 result: ", "Chromosome\tStart\tEnd\tFeature\tMR[sample1]\n" +
                "Case3\t4\t21\tannotation4\t52.7778\n" +
                "Case3\t11\t26\tannotation6\t35.8209\n" +
                "Case3\t4\t30\tannotation3\t43.0000\n" +
                "Case3\t11\t30\tannotation5\t28.1250\n", stringWriter.getBuffer().toString());
    }

    @Test
    public void testCase4() {
        String[] samples = new String[]{"sample1", "sample2", "sample3"};
        int[] positions = new int[]{6, 8, 14, 16};
        int[][] C = {{5, 3, 9, 8}, {4, 6, 3, 2}, {8, 3, 8, 9}};
        int[][] Cm = {{9, 7, 1, 5}, {9, 7, 9, 3}, {2, 3, 2, 8}};
        testSupport = new MethylCountProviderTestSupport(samples, positions, "Case2", C, Cm);
        final StringWriter stringWriter = new StringWriter();
        VCFAveragingWriter testWriter = new VCFAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 4 result: ",
                "Chromosome\tStart\tEnd\tFeature\tMR[sample1]\tMR[sample2]\tMR[sample3]\n" +
                        "Case2\t5\t9\tannotation1\t66.6667\t61.5385\t31.2500\n" +
                        "Case2\t13\t17\tannotation2\t26.0870\t70.5882\t37.0370\n", stringWriter.getBuffer().toString());
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
        VCFAveragingWriter testWriter = new VCFAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        int[] a = {0, 0, 0};
        testWriter.setSampleIndexToGroupIndex(a);
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.writeRecord();
        testWriter.close();
        assertEquals("Test Case 5 result: ",
                "Chromosome\tStart\tEnd\tFeature\tMR[sample1]\tMR[sample2]\tMR[sample3]\tMR[group1]\n" +
                        "Case2\t5\t9\tannotation1\t66.6667\t61.5385\t31.2500\t56.0606\n" +
                        "Case2\t13\t17\tannotation2\t26.0870\t70.5882\t37.0370\t41.7910\n", stringWriter.getBuffer().toString());

    }

    @Test
    public void testCase6() {
        String[] groups = new String[]{"group1", "group2"};
        String[] samples = new String[]{"sample1", "sample2", "sample3", "sample4", "sample5", "sample6"};
        int[] positions = new int[]{5, 7, 12, 14, 19, 24, 26, 29};
        int[][] C = {{5,3,4,8,5,2,8,7}, {5,3,9,8,9,8,4,9}, {5,3,9,8,7,9,8,0}, {5,3,9,8,9,9,8,0},
                {5,3,9,8,2,9,8,0},{0,3,5,8,8,7,8,0}};
        int[][] Cm = {{9,7,1,5,3,7,0,3},{9,3,1,7,4,7,0,3},{9,6,1,2,8,7,8,3},{9,1,1,9,2,7,8,3},
                {0,7,1,8,8,7,8,3},{8,7,1,5,3,7,8,3}};
        testSupport = new MethylCountProviderTestSupport(groups, samples, positions, "Case3", C, Cm);

        final StringWriter stringWriter = new StringWriter();
        VCFAveragingWriter testWriter = new VCFAveragingWriter(stringWriter, genome, testSupport);
        testWriter.setAnnotationFilename("test-data/vcf-averaging/annotations-1.tsv");
        int[] a={0,0,0,1,1,1};
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
        assertEquals("Test Case 6 result: ", "Chromosome\tStart\tEnd\tFeature\tMR[sample1]\tMR[sample2]" +
                "\tMR[sample3]\tMR[sample4]\tMR[sample5]\tMR[sample6]\tMR[group1]\tMR[group2]\n" +
                "Case3\t4\t21\tannotation4\t59.3750\t48.4848\t60.5263\t41.3793\t60.0000\t62.0690\t56.3107\t54.2169\n" +
                "Case3\t11\t26\tannotation6\t37.2093\t33.3333\t38.8060\t38.5714\t47.0588\t40.0000\t36.5269\t41.9192\n" +
                "Case3\t4\t30\tannotation3\t45.4545\t38.2022\t47.3118\t43.9560\t48.8372\t51.8519\t43.6293\t48.0620\n" +
                "Case3\t11\t30\tannotation5\t32.1429\t29.7297\t26.0870\t43.3333\t41.3793\t40.9091\t29.5455\t41.9753\n", stringWriter.getBuffer().toString());

    }
}
