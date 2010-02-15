/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: xutao
 * Date: May 17, 2009
 * Time: 1:02:50 AM
 * To change this template use File | Settings | File Templates.
 */
public class TestAnnotationCount {
    private AnnotationCount annotationCount;
    private static final Log LOG = LogFactory.getLog(TestAnnotationCount.class);

    @Before
    public void setUp() throws IOException {
        annotationCount = new AnnotationCount();
    }

    @Test
    public void testAnnotCountDepth() {
        annotationCount.baseCounter.startPopulating();
        annotationCount.populate(3, 8);
        annotationCount.populate(9, 10);
        annotationCount.populate(5, 12);
        annotationCount.populate(3, 7);
        annotationCount.populate(8, 12);

        final int[] trueCount = {0, 0, 0, 2, 2, 3, 3, 3, 3, 3, 3, 2, 2, 0, 0, 0};
        annotationCount.baseCounter.accumulate();
        annotationCount.baseCounter.baseCount(); // final algorithm for base count without writer
        System.out.println("count perbase" + annotationCount.baseCounter.countPerBase);
        System.out.println("count keys " + annotationCount.baseCounter.countKeys);
        for (int i = 0; i <= 15; i++) {
            for (int j = 0; j <= 15; j++) {
                float sum = 0;
                for (int p = i; p <= j; p++) {
                    sum += trueCount[p];
                }
                final float expDepth = sum == 0 ? 0 : sum / (j - i + 1);
                assertEquals(expDepth, annotationCount.averageReadsPerPosition(i, j), 0.01);
            }
        }
    }

    @Test
    public void testAnnotCountOverlapCount() {
        annotationCount.baseCounter.startPopulating();
        annotationCount.populate(3, 8);
        annotationCount.populate(9, 10);
        annotationCount.populate(5, 12);
        annotationCount.populate(3, 7);
        annotationCount.populate(8, 12);
        annotationCount.populate(15, 18);
        annotationCount.sortReads();
        annotationCount.baseCounter.accumulate();
        annotationCount.baseCounter.baseCount(); // final algorithm for base count without writer
        for(int i=-1; i<=30; i++) {
            System.out.println(i + " " + annotationCount.countReadsPartiallyOverlappingWithInterval(i, i));
        }
        System.out.println("overlapping count 3, 4 " + annotationCount.countReadsPartiallyOverlappingWithInterval(3, 4));
        System.out.println("overlapping count 2, 3 " + annotationCount.countReadsPartiallyOverlappingWithInterval(2, 3));
        System.out.println("overlapping count 3, 8 " + annotationCount.countReadsPartiallyOverlappingWithInterval(3, 8));
        System.out.println("overlapping count 11, 12 " + annotationCount.countReadsPartiallyOverlappingWithInterval(11, 12));
        System.out.println("overlapping count 9, 10 " + annotationCount.countReadsPartiallyOverlappingWithInterval(9, 10));
        System.out.println("overlapping count 8, 9 " + annotationCount.countReadsPartiallyOverlappingWithInterval(8, 9));
        System.out.println("overlapping count 8, 8 " + annotationCount.countReadsPartiallyOverlappingWithInterval(8, 8));
        System.out.println("overlapping count 5, 6 " + annotationCount.countReadsPartiallyOverlappingWithInterval(5, 6));
        System.out.println("overlapping count 3, 3 " + annotationCount.countReadsPartiallyOverlappingWithInterval(3, 3));
        System.out.println("overlapping count 3, 12 " + annotationCount.countReadsPartiallyOverlappingWithInterval(3, 12));
        System.out.println("overlapping count -1, 2 " + annotationCount.countReadsPartiallyOverlappingWithInterval(-1, 2));
        System.out.println("overlapping count 0, 45 " + annotationCount.countReadsPartiallyOverlappingWithInterval(0, 45));
        System.out.println("overlapping count 13, 15 " + annotationCount.countReadsPartiallyOverlappingWithInterval(13, 15));

    }

    @Test
    public void testAnnotCountInsideCount() {
        annotationCount.baseCounter.startPopulating();
        annotationCount.populate(3, 8 );
        annotationCount.populate(9, 10 );
        annotationCount.populate(5, 12 );
        annotationCount.populate(3, 7 );
        annotationCount.populate(8, 12 );
        annotationCount.populate(15, 18 );
        annotationCount.sortReads();
        annotationCount.baseCounter.accumulate();
        annotationCount.baseCounter.baseCount(); // final algorithm for base count without writer
        for(int i=-1; i<=30; i++) {
            System.out.println(i + " " + annotationCount.countReadsStriclyWithinInterval(i, i));
        }
        System.out.println("overlapping count 3, 4 " + annotationCount.countReadsStriclyWithinInterval(3, 4));
        System.out.println("overlapping count 15, 17 " + annotationCount.countReadsStriclyWithinInterval(15, 17));
        System.out.println("overlapping count 9, 18 " + annotationCount.countReadsStriclyWithinInterval(9, 18));
        System.out.println("overlapping count 3, 8 " + annotationCount.countReadsStriclyWithinInterval(3, 8));
        System.out.println("overlapping count 11, 12 " + annotationCount.countReadsStriclyWithinInterval(11, 12));
        System.out.println("overlapping count 9, 10 " + annotationCount.countReadsStriclyWithinInterval(9, 10));
        System.out.println("overlapping count 8, 9 " + annotationCount.countReadsStriclyWithinInterval(8, 9));
        System.out.println("overlapping count 8, 8 " + annotationCount.countReadsStriclyWithinInterval(8, 8));
        System.out.println("overlapping count 5, 6 " + annotationCount.countReadsStriclyWithinInterval(5, 6));
        System.out.println("overlapping count 3, 3 " + annotationCount.countReadsStriclyWithinInterval(3, 3));
        System.out.println("overlapping count 3, 12 " + annotationCount.countReadsStriclyWithinInterval(3, 12));
        System.out.println("overlapping count -1, 2 " + annotationCount.countReadsStriclyWithinInterval(-1, 2));
        System.out.println("overlapping count 0, 45 " + annotationCount.countReadsStriclyWithinInterval(0, 45));
        System.out.println("overlapping count 13, 15 " + annotationCount.countReadsStriclyWithinInterval(13, 15));
    }

    @Test
    public void testGeneExpression() {
        annotationCount.baseCounter.startPopulating();
        annotationCount.populate(3, 8 );
        annotationCount.populate(9, 10 );
        annotationCount.populate(5, 12 );
        annotationCount.populate(3, 7 );
        annotationCount.populate(8, 12 );
        annotationCount.populate(15, 18 );
        annotationCount.sortReads();
        annotationCount.baseCounter.accumulate();
        annotationCount.baseCounter.baseCount(); // final algorithm for base count without writer
    }

    // TODO @Test
    public void testReadingAnnotation() throws IOException {
//        CompactAlignmentToAnnotationCountsMode b =new CompactAlignmentToAnnotationCountsMode();
//        Object2ObjectMap<String, ObjectList<Annotation>> x= b.readAnnotations("C:\\cygwin\\My Dropbox\\icb-hadoop\\data\\human_exon_annotation_biomart_NCBI36.txt");
//        int sum=0;
//        for(String chro: x.keySet()){
//            ObjectList <Annotation> annots= x.get(chro);
//            System.out.println(chro+"   "+annots.size());
//            for(Annotation annot: annots){
//                System.out.println("id  "+annot.id);
//                sum+=annot.segments.size();
//            }
//        }
//            System.out.println("sum "+sum);
    }
}
