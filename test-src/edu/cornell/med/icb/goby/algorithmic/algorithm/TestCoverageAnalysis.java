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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.counts.AnyTransitionCountsIterator;
import edu.cornell.med.icb.goby.counts.CountsReaderI;
import edu.cornell.med.icb.goby.counts.CountsReaderTestSupport;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.lang.MutableString;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: May 21, 2011
 *         Time: 11:36:33 AM
 */
public class TestCoverageAnalysis {
    @Test
    public void constantDepth() throws IOException {
        CoverageAnalysis analysis = new CoverageAnalysis();
        final int[] annotLengths = {10};
        int[] annotCounts = {1};
        CountsReaderI annotations = new CountsReaderTestSupport(annotLengths, annotCounts);

        int[] lengths = {10};
        int[] counts = {4};
        CountsReaderI reader = new CountsReaderTestSupport(lengths, counts);

        analysis.process(annotations,reader);

        assertEquals(4d, analysis.getAverageDepth(), .01);
    }

    @Test
    public void constantDepth2() throws IOException {
        CoverageAnalysis analysis = new CoverageAnalysis();
        int[] annotLengths = {20, 10, 0};
        int[] annotCounts = {0, 1, 0};
        CountsReaderI annotations = new CountsReaderTestSupport(annotLengths, annotCounts);

        int[] lengths = {20, 10};
        int[] counts = {0, 4};
        CountsReaderI reader = new CountsReaderTestSupport(lengths, counts);

        analysis.process(annotations, reader);

        assertEquals(4d, analysis.getAverageDepth(), .01);
        assertEquals(4d, analysis.getAnnotationAverageDepth(), .01);
    }

    @Test
    public void annotationDepth() throws IOException {
        CoverageAnalysis analysis = new CoverageAnalysis();
        int[] annotLengths = {2, 3, 2};
        int[] annotCounts = {0, 1, 0};
        CountsReaderI annotations = new CountsReaderTestSupport(annotLengths, annotCounts);

        int[] lengths = {2, 3, 2};
        int[] counts = {2, 4, 2};
        CountsReaderI reader = new CountsReaderTestSupport(lengths, counts);

        analysis.process(annotations, reader);

        assertEquals((2.0d * 4 + 4 * 3) / 7.0d, analysis.getAverageDepth(), 0.01d);
        assertEquals(4.0d, analysis.getAnnotationAverageDepth(), 0.01);
        assertEquals(2.0d, analysis.getNotAnnotationAverageDepth(), 0.01);
        analysis.estimateStatistics();
        assertEquals(0.0d, analysis.getNumBasesWithDepthAtLeast(5), 0.1);
        assertEquals(3 * 4.0d, analysis.getNumBasesWithDepthAtLeast(4), 0.1);
        assertEquals(3.0d, analysis.getNumSitesWithDepthAtLeast(4), 0.1);

    }

    public void testIterator() throws IOException {


        // (length, count) (2,0) (8,1) (1,0)
        int[] annotLengths = {2, 8, 1};
        int[] annotCounts = {0, 1, 0};
        CountsReaderI annotations = new CountsReaderTestSupport(annotLengths, annotCounts);

        // (length, count) (2,2) (3,7) (3,4) (3,2)
        int[] lengths = {2, 3, 3, 3};
        int[] counts = {2, 7, 4, 2};
        CountsReaderI reader = new CountsReaderTestSupport(lengths, counts);
        AnyTransitionCountsIterator orIterator = new AnyTransitionCountsIterator(reader, annotations);

        orIterator.nextTransition();
        assertEquals(0, orIterator.getCount(1));
        assertEquals(2, orIterator.getCount(0));
        assertEquals(2, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(1, orIterator.getCount(1));
        assertEquals(7, orIterator.getCount(0));
        assertEquals(3, orIterator.getPosition());
        orIterator.nextTransition();
        assertEquals(1, orIterator.getCount(1));
        assertEquals(7, orIterator.getCount(0));
        assertEquals(5, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(1, orIterator.getCount(1));
        assertEquals(4, orIterator.getCount(0));
        assertEquals(8, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(1, orIterator.getCount(1));
        assertEquals(2, orIterator.getCount(0));
        assertEquals(10, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(0, orIterator.getCount(1));
        assertEquals(2, orIterator.getCount(0));
        assertEquals(11, orIterator.getPosition());
        assertEquals(1, orIterator.getLength());

    }

    @Test
    public void capturedDepthAtPercentile() throws IOException {
        CoverageAnalysis analysis = new CoverageAnalysis();
        // (length, count) (2,0) (8,1) (1,0)
        // TODO fix problem when annotations do not end at the same position as counts.
        int[] annotLengths = {2, 8, 1};
        int[] annotCounts = {0, 1, 0};
        CountsReaderI annotations = new CountsReaderTestSupport(annotLengths, annotCounts);

        // (length, count) (2,2) (3,7) (3,4) (3,2)
        int[] lengths = {2, 3, 3, 3};
        int[] counts = {2, 7, 4, 2};
        CountsReaderI reader = new CountsReaderTestSupport(lengths, counts);

        analysis.process(annotations, reader);
        LongArrayList inside = analysis.getDepthTallySitesInAnnotation();
        assertEquals(3l, (long) inside.get(7));
        assertEquals(0l, (long) inside.get(5));
        assertEquals(3l, (long) inside.get(4));
        assertEquals(0l, (long) inside.get(3));
        assertEquals(2l, (long) inside.get(2));
        analysis.estimateStatistics();

        long[] totalCumul = analysis.getCumulativeSitesTotal();

        assertEquals(11l, (long) totalCumul[0]);
        assertEquals(11l, (long) totalCumul[1]);
        assertEquals(11l, (long) totalCumul[2]);
        assertEquals(6l, (long) totalCumul[3]);
        assertEquals(6l, (long) totalCumul[4]);
        assertEquals(3l, (long) totalCumul[5]);
        assertEquals(3l, (long) totalCumul[6]);
        assertEquals(3l, (long) totalCumul[7]);


        long[] insideCumul = analysis.getCumulativeSitesCaptured();

        assertEquals(8l, (long) insideCumul[0]);
        assertEquals(8l, (long) insideCumul[1]);
        assertEquals(8l, (long) insideCumul[2]);
        assertEquals(6l, (long) insideCumul[3]);
        assertEquals(6l, (long) insideCumul[4]);
        assertEquals(3l, (long) insideCumul[5]);
        assertEquals(3l, (long) insideCumul[6]);
        assertEquals(3l, (long) insideCumul[7]);


        long[] outsideCumul = analysis.getCumulativeSitesNotCaptured();

        assertEquals(3l, (long) outsideCumul[0]);
        assertEquals(3l, (long) outsideCumul[1]);
        assertEquals(3l, (long) outsideCumul[2]);
        assertEquals(0l, (long) outsideCumul[3]);
        assertEquals(0l, (long) outsideCumul[4]);
        assertEquals(0l, (long) outsideCumul[5]);
        assertEquals(0l, (long) outsideCumul[6]);
        assertEquals(0l, (long) outsideCumul[7]);


        assertEquals(8d, analysis.getNumSitesCapturedWithDepthAtLeast(0), .1);
        assertEquals(8d, analysis.getNumSitesCapturedWithDepthAtLeast(1), .1);
        assertEquals(6d, analysis.getNumSitesCapturedWithDepthAtLeast(4), .1);
        assertEquals(3d, analysis.getNumSitesCapturedWithDepthAtLeast(7), .1);
        assertEquals(0d, analysis.getNumSitesCapturedWithDepthAtLeast(8), .1);

        assertEquals(11d, analysis.getNumSitesWithDepthAtLeast(0), .1);
        assertEquals(11d, analysis.getNumSitesWithDepthAtLeast(1), .1);
        assertEquals(11d, analysis.getNumSitesWithDepthAtLeast(2), .1);
        assertEquals(6d, analysis.getNumSitesWithDepthAtLeast(3), .1);
        assertEquals(6d, analysis.getNumSitesWithDepthAtLeast(4), .1);
        assertEquals(3d, analysis.getNumSitesWithDepthAtLeast(7), .1);

        assertEquals(5d, analysis.depthAtPercentile(.5), .1);
        assertEquals(3d, analysis.depthCapturedAtPercentile(.9), .1);
        assertEquals(7d, analysis.depthCapturedAtPercentile(.1), .1);

    }


    //TODO this test must pass without the final (1,0) in annotation.
  @Test
    public void problemWithAnnotationsTooShort() throws IOException {
        CoverageAnalysis analysis = new CoverageAnalysis();
        // (length, count) (2,0) (8,1) (1,0)

        int[] annotLengths = {2, 8};
        int[] annotCounts = {0, 1};
        CountsReaderI annotations = new CountsReaderTestSupport(annotLengths, annotCounts);

        // (length, count) (2,2) (3,7) (3,4) (3,2)
        int[] lengths = {2, 3, 3, 3};
        int[] counts = {2, 7, 4, 2};
        CountsReaderI reader = new CountsReaderTestSupport(lengths, counts);

        analysis.process(annotations, reader);
        LongArrayList inside = analysis.getDepthTallySitesInAnnotation();
        assertEquals(3l, (long) inside.get(7));
        assertEquals(0l, (long) inside.get(5));
        assertEquals(3l, (long) inside.get(4));
        assertEquals(0l, (long) inside.get(3));
        assertEquals(2l, (long) inside.get(2));
        analysis.estimateStatistics();

        long[] totalCumul = analysis.getCumulativeSitesTotal();

        assertEquals(11l, (long) totalCumul[0]);
        assertEquals(11l, (long) totalCumul[1]);
        assertEquals(11l, (long) totalCumul[2]);
        assertEquals(6l, (long) totalCumul[3]);
        assertEquals(6l, (long) totalCumul[4]);
        assertEquals(3l, (long) totalCumul[5]);
        assertEquals(3l, (long) totalCumul[6]);
        assertEquals(3l, (long) totalCumul[7]);


        long[] insideCumul = analysis.getCumulativeSitesCaptured();

        assertEquals(8l, (long) insideCumul[0]);
        assertEquals(8l, (long) insideCumul[1]);
        assertEquals(8l, (long) insideCumul[2]);
        assertEquals(6l, (long) insideCumul[3]);
        assertEquals(6l, (long) insideCumul[4]);
        assertEquals(3l, (long) insideCumul[5]);
        assertEquals(3l, (long) insideCumul[6]);
        assertEquals(3l, (long) insideCumul[7]);


        long[] outsideCumul = analysis.getCumulativeSitesNotCaptured();

        assertEquals(3l, (long) outsideCumul[0]);
        assertEquals(3l, (long) outsideCumul[1]);
        assertEquals(3l, (long) outsideCumul[2]);
        assertEquals(0l, (long) outsideCumul[3]);
        assertEquals(0l, (long) outsideCumul[4]);
        assertEquals(0l, (long) outsideCumul[5]);
        assertEquals(0l, (long) outsideCumul[6]);
        assertEquals(0l, (long) outsideCumul[7]);
    }


    // test may be incorrect
    public void iteratorProblemAnnotationTooShort() throws IOException {


        // (length, count) (2,0) (8,1) (1,0)
        int[] annotLengths = {2, 8};
        int[] annotCounts = {0, 1};
        CountsReaderI annotations = new CountsReaderTestSupport(annotLengths, annotCounts);

        // (length, count) (2,2) (3,7) (3,4) (3,2)
        int[] lengths = {2, 3, 3, 3};
        int[] counts = {2, 7, 4, 2};
        CountsReaderI reader = new CountsReaderTestSupport(lengths, counts);
        AnyTransitionCountsIterator orIterator = new AnyTransitionCountsIterator(reader, annotations);

        orIterator.nextTransition();
        assertEquals(0, orIterator.getCount(1));
        assertEquals(2, orIterator.getCount(0));
        assertEquals(2, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(1, orIterator.getCount(1));
        assertEquals(7, orIterator.getCount(0));
        assertEquals(5, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(1, orIterator.getCount(1));
        assertEquals(4, orIterator.getCount(0));
        assertEquals(8, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(1, orIterator.getCount(1));
        assertEquals(2, orIterator.getCount(0));
        assertEquals(10, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(0, orIterator.getCount(1));
        assertEquals(2, orIterator.getCount(0));
        assertEquals(11, orIterator.getPosition());
        assertEquals(1, orIterator.getLength());

    }

    // test may be incorrect
    public void iteratorLargeAnnotation() throws IOException {


        // (length, count) (2,0) (8,1) (1,0)
        int[] annotLengths = {2, 5};
        int[] annotCounts = {0, 1};
        CountsReaderI annotations = new CountsReaderTestSupport(annotLengths, annotCounts);

        // (length, count) (2,2) (3,7) (3,4) (3,2)
        int[] lengths = {2, 3, 3, 3};
        int[] counts = {2, 7, 4, 2};
        CountsReaderI reader = new CountsReaderTestSupport(lengths, counts);
        AnyTransitionCountsIterator orIterator = new AnyTransitionCountsIterator(reader, annotations);

        orIterator.nextTransition();
        assertEquals(0, orIterator.getCount(1));
        assertEquals(2, orIterator.getCount(0));
        assertEquals(2, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(1, orIterator.getCount(1));
        assertEquals(7, orIterator.getCount(0));
        assertEquals(5, orIterator.getPosition());

        orIterator.nextTransition();
        assertEquals(1, orIterator.getCount(1));
        assertEquals(4, orIterator.getCount(0));
        assertEquals(7, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(0, orIterator.getCount(1));
        assertEquals(4, orIterator.getCount(0));
        assertEquals(8, orIterator.getPosition());

        orIterator.nextTransition();

        assertEquals(0, orIterator.getCount(1));
        assertEquals(2, orIterator.getCount(0));
        //    assertEquals(11, orIterator.getPosition());
        assertEquals(3, orIterator.getLength());
    }

    @Test
    public void small() throws IOException {


        // (length, count) (2,0) (8,1) (1,0)


        assertEquals("(1,0)(3,1)", outputLengthCount("(1,0)(3,1)"));
        assertEquals("(1,0)(4,1)", outputPositionCount("(1,0)(3,1)"));
    }

    
    @Test
    public void twoFlats() throws IOException {


        // (length, count) (2,0) (8,1) (1,0)


        String expected = "(2,0)(1,1)(2,0)(2,1)(2,2)(1,1)(1,2)(2,1)(1,2)(1,1)";

        String produced =
                outputLengthCount("(2,0)(1,1)(4,0)(2,1)(1,0)(1,1)(2,0)(1,1)(3,0)(1,1)",
                "(5,0)(10,1)");
        //    assertEquals(expected, produced);


        assertEquals("(2,0)(3,1)(5,0)(7,1)(9,2)(10,1)(11,2)(13,1)(14,2)(15,1)(17,0)(18,1)",

                outputPositionCount("(2,0)(1,1)(4,0)(2,1)(1,0)(1,1)(2,0)(1,1)(3,0)(1,1)",
                "(5,0)(10,1)"));
        /*
        int[][] expected = {
                // {0, 0, 0},
                {1, 0, 0},
                {2, 1, 0},
                {3, 0, 0},
                //          {4, 0, 0},
                {5, 0, 1},
                //        {6, 0, 1},
                {7, 1, 1},
                //      {8, 1, 1},
                {9, 0, 1},
                {10, 1, 1},
                {11, 0, 1},
                //    {12, 0, 1},
                {13, 1, 1},
                {14, 0, 1},
                {15, 0, 0},
                //  {16, 0, 0},
                {17, 1, 0}};


        orIterator = new AnyTransitionCountsIterator(reader0, reader1);

        int index = 0;
        while (orIterator.hasNextTransition()) {
            orIterator.nextTransition();
            int position = orIterator.getPosition();
            System.out.printf("index=%d position=%d %n", index, position);
            assertEquals(expected[index][0], position);
            for (int readerIndex = 0; readerIndex < 2; readerIndex++) {
                assertEquals(expected[index][readerIndex + 1], orIterator.getCount(readerIndex));
            }
            index++;
        }
        */
    }
      @Test
    public void testNegativeStart() throws IOException {

          assertEquals("(2,0)(4,1)(5,2)(62,1)",
                outputPositionCount("(4,0)(1,1)",
                "(2,0)(60,1)"));
    }
       @Test
    public void testAnnotation() throws IOException {

          assertEquals("(32271285,1)(32271290,1)(32271295,1)",
                outputPositionCount("(32271280,0)(770290,1)",
                "(32271285,0)(1,1)(5,0)(1,1)(5,0)(1,1)"));
    }


    private String outputLengthCount(final String... formats) throws IOException {
        CountsReaderI readers[] = new CountsReaderI[formats.length];
        int i = 0;
        for (String format : formats) {
            readers[i] = new CountsReaderTestSupport(format);
            i++;
        }
        AnyTransitionCountsIterator orIterator;
        MutableString result = new MutableString();
        orIterator = new AnyTransitionCountsIterator(readers);
        while (orIterator.hasNextTransition()) {
            orIterator.nextTransition();

            result.append(String.format("(%d,%d)", orIterator.getLength(), orIterator.getCount()));

        }
        return result.toString();
    }

    private String outputPositionCount(final String... formats) throws IOException {
        CountsReaderI readers[] = new CountsReaderI[formats.length];
        int i = 0;
        for (String format : formats) {
            readers[i] = new CountsReaderTestSupport(format);
            i++;
        }
        AnyTransitionCountsIterator orIterator;
        MutableString result = new MutableString();
        orIterator = new AnyTransitionCountsIterator(readers);
        while (orIterator.hasNextTransition()) {
            orIterator.nextTransition();

            result.append(String.format("(%d,%d)", orIterator.getPosition(), orIterator.getCount()));

        }
        return result.toString();
    }


    //TODO reactivate these tests
    public void fourFlats() throws IOException {


        // (length, count) (2,0) (8,1) (1,0)

        CountsReaderI reader0 = new CountsReaderTestSupport("(5,0)(1,1)(3,0)(1,1)(3,0)(1,1)(1,0)(1,1)(1,0)(1,1) ");
        CountsReaderI reader1 = new CountsReaderTestSupport("(2,0)(6,1)");
        CountsReaderI reader2 = new CountsReaderTestSupport("(3,0)(7,1) ");
        CountsReaderI reader3 = new CountsReaderTestSupport("(5,0)(10,1) ");
        AnyTransitionCountsIterator orIterator;

        /* orIterator = new AnyTransitionCountsIterator(reader0, reader1, reader2, reader3);
       while (orIterator.hasNextTransition()) {
           orIterator.nextTransition();

           System.out.printf("position=%d count=%d length=%d%n", orIterator.getPosition(), orIterator.getCount(), orIterator.getLength());

       }
        */
        int[][] expected = {
                //  {0, 0, 0, 0, 0},
                //  {1, 0, 0, 0, 0},
                {2, 0, 1, 0, 0},
                {3, 0, 1, 1, 0},
                {4, 0, 1, 1, 0},
                {5, 1, 1, 1, 1},
                {6, 0, 1, 1, 1},
                {7, 0, 1, 1, 1},
                {8, 0, 0, 1, 1},
                {9, 1, 0, 1, 1},
                {10, 0, 0, 0, 1},
                {11, 0, 0, 0, 1},
                {12, 0, 0, 0, 1},
                {13, 1, 0, 0, 1},
                {14, 0, 0, 0, 1},
                {15, 1, 0, 0, 0},
                {16, 0, 0, 0, 0},
                {17, 1, 0, 0, 0}};


        orIterator = new AnyTransitionCountsIterator(reader0, reader1, reader2, reader3);

        int index = 0;
        while (orIterator.hasNextTransition()) {
            orIterator.nextTransition();
            int position = orIterator.getPosition();
            System.out.printf("index=%d position=%d %n", index, position);
            assertEquals(expected[index][0], position);
            for (int readerIndex = 0; readerIndex < 4; readerIndex++) {
                assertEquals(expected[index][readerIndex + 1], orIterator.getCount(readerIndex));
            }
            index++;
        }

    }
}


