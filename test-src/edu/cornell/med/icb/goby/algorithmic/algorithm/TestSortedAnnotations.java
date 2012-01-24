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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.Segment;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.junit.*;

import java.io.IOException;
import java.io.StringReader;
import java.util.concurrent.atomic.AtomicReference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * @author Nyasha Chambwe
 *         Date: 1/10/12
 *         Time: 2:31 PM
 */

public class TestSortedAnnotations {

    SortedAnnotations testObj = new SortedAnnotations();

    Annotation[] annotations;

    String[] sequences = {
            "Case1",
            "Case2",
            "Case3",
            "Chr1",
            "Chr2"
    };

    RandomAccessSequenceTestSupport genome = new RandomAccessSequenceTestSupport(sequences){
        @Override
        public int getReferenceIndex(String referenceId) {

           if(referenceId.equals("Case1")){return 0;}
           if(referenceId.equals("Case2")){return 1;}
           if(referenceId.equals("Case3")){return 2;}
           if(referenceId.equals("Chr1")){return 3;}
           if(referenceId.equals("Chr2")){return 4;}

            return  -1;
        }

        @Override
        public String getReferenceName(int index) {
           return sequences[index];
        }
    };

   @Before
   public void setUp() throws IOException {
        
        StringReader reader=new StringReader("Chromosome\tStrand\tTranscriptID\tSegmentID\tStart\tEnd\n"+
                "Case1\t+\tannotation0\t1\t3\t11\n" +
                "Case1\t+\tannotation1\t2\t6\t10\n" +
                "Case1\t+\tannotation2\t3\t7\t13\n" +
                "Case1\t+\tannotation3\t4\t2\t4\n" +
                "Case1\t+\tannotation3\t5\t6\t8\n" +
                "Case1\t+\tannotation3\t6\t11\t14\n" +
                "Case2\t+\tannotation6\t7\t5\t8\n" +
                "Case2\t+\tannotation7\t8\t16\t19\n" +
                "Case3\t+\tannotation8\t9\t2\t5\n" +
                "Case3\t+\tannotation8\t10\t13\t17\n" +
                "Case3\t+\tannotation8\t11\t26\t29\n" +
                "Case3\t+\tannotation11\t12\t6\t9\n" +
                "Case3\t+\tannotation11\t13\t17\t20\n" +
                "Case3\t+\tannotation11\t14\t28\t29\n"+
                "Chr1\t+\tannotation12\t15\t4\t7\n"+
                "Chr2\t+\tannotation13\t15\t4\t7\n"
        );
        testObj.setGenome(genome);
        testObj.loadAnnotations(reader);

        System.out.println("Annotations loaded");
    }

    public ObjectArrayList<Annotation> getCurrentAnnotations(int refIndex, int pos){
        testObj.hasOverlappingAnnotations(refIndex, pos);
        return testObj.currentAnnotations();
    }

    @Test
    public void TestSortedAnnotationsCase1() {

        ObjectArrayList<Annotation> currentAnnotations;
        int refIndex= genome.getReferenceIndex("Case1");
        
        assertTrue(testObj.currentAnnotations().isEmpty());

        currentAnnotations= getCurrentAnnotations(refIndex, 0);
        assertEquals("Position 0: ", "[]", testObj.currentAnnotations().toString());

        // return set and make sure it's empty
        assertTrue(testObj.currentAnnotations().isEmpty());


        currentAnnotations= getCurrentAnnotations(refIndex, 1);
        assertEquals("Position 1:", "[]", currentAnnotations.toString());

        assertTrue(testObj.currentAnnotations().isEmpty());

        currentAnnotations= getCurrentAnnotations(refIndex, 2);
        assertEquals("Position 2:", "[[ Case1 annotation3 2-4 6-8 11-14  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 3);
        assertEquals("Position 3:", "[[ Case1 annotation3 2-4 6-8 11-14  ], [ Case1 annotation0 3-11  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 4);
        assertEquals("Position 4:", "[[ Case1 annotation3 2-4 6-8 11-14  ], [ Case1 annotation0 3-11  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 5);
        assertEquals("Position 5:", "[[ Case1 annotation0 3-11  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 6);
        assertEquals("Position 6:", "[[ Case1 annotation3 2-4 6-8 11-14  ], [ Case1 annotation0 3-11  ], [ Case1 annotation1 6-10  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 7);
        assertEquals("Position 7:", "[[ Case1 annotation3 2-4 6-8 11-14  ], [ Case1 annotation0 3-11  ], [ Case1 annotation1 6-10  ], [ Case1 annotation2 7-13  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 8);
        assertEquals("Position 8:", "[[ Case1 annotation3 2-4 6-8 11-14  ], [ Case1 annotation0 3-11  ], [ Case1 annotation1 6-10  ], [ Case1 annotation2 7-13  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 9);
        assertEquals("Position 9:", "[[ Case1 annotation0 3-11  ], [ Case1 annotation1 6-10  ], [ Case1 annotation2 7-13  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 10);
        assertEquals("Position 10:", "[[ Case1 annotation0 3-11  ], [ Case1 annotation1 6-10  ], [ Case1 annotation2 7-13  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 11);
        assertEquals("Position 11:", "[[ Case1 annotation3 2-4 6-8 11-14  ], [ Case1 annotation0 3-11  ], [ Case1 annotation2 7-13  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 12);
        assertEquals("Position 12:", "[[ Case1 annotation3 2-4 6-8 11-14  ], [ Case1 annotation2 7-13  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 13);
        assertEquals("Position 13:", "[[ Case1 annotation3 2-4 6-8 11-14  ], [ Case1 annotation2 7-13  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 14);
        assertEquals("Position 14:", "[[ Case1 annotation3 2-4 6-8 11-14  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 15);
        assertEquals("Position 15:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 0);
        assertEquals("Position 0:", "[]", currentAnnotations.toString());
        assertTrue(currentAnnotations.isEmpty());

       
    }

    @Test
    public void TestSortedAnnotationsCase2() {

        int refIndex= genome.getReferenceIndex("Case2");
        ObjectArrayList<Annotation> currentAnnotations;

        currentAnnotations= getCurrentAnnotations(refIndex, 0);
        assertEquals("Position 0:", "[]", currentAnnotations.toString());
        assertTrue(currentAnnotations.isEmpty());

        currentAnnotations= getCurrentAnnotations(refIndex, 1);
        assertEquals("Position 1:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 2);
        assertEquals("Position 2:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 3);
        assertEquals("Position 3:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 4);
        assertEquals("Position 4:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 5);
        assertEquals("Position 5:", "[[ Case2 annotation6 5-8  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 6);
        assertEquals("Position 6:", "[[ Case2 annotation6 5-8  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 7);
        assertEquals("Position 7:", "[[ Case2 annotation6 5-8  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 8);
        assertEquals("Position 8:", "[[ Case2 annotation6 5-8  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 9);
        assertEquals("Position 9:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 10);
        assertEquals("Position 10:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 11);
        assertEquals("Position 11:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 12);
        assertEquals("Position 12:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 13);
        assertEquals("Position 13:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 14);
        assertEquals("Position 14:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 15);
        assertEquals("Position 15:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 16);
        assertEquals("Position 16:", "[[ Case2 annotation7 16-19  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 17);
        assertEquals("Position 17:", "[[ Case2 annotation7 16-19  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 18);
        assertEquals("Position 18:", "[[ Case2 annotation7 16-19  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 19);
        assertEquals("Position 19:", "[[ Case2 annotation7 16-19  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 20);
        assertEquals("Position 20:", "[]", currentAnnotations.toString());

    }

    @Test
    public void TestSortedAnnotationsCase3() {

            int refIndex= genome.getReferenceIndex("Case3");
            ObjectArrayList<Annotation> currentAnnotations;

            currentAnnotations= getCurrentAnnotations(refIndex, 0);
            assertEquals("Position 0:", "[]", currentAnnotations.toString());
            assertTrue(currentAnnotations.isEmpty());

            currentAnnotations= getCurrentAnnotations(refIndex, 1);
            assertEquals("Position 1:", "[]", currentAnnotations.toString());


         currentAnnotations= getCurrentAnnotations(refIndex, 2);
        assertEquals("Position 2:", "[[ Case3 annotation8 2-5 13-17 26-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 3);
        assertEquals("Position 3:", "[[ Case3 annotation8 2-5 13-17 26-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 4);
        assertEquals("Position 4:", "[[ Case3 annotation8 2-5 13-17 26-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 5);
        assertEquals("Position 5:", "[[ Case3 annotation8 2-5 13-17 26-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 6);
        assertEquals("Position 6:", "[[ Case3 annotation11 6-9 17-20 28-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 7);
        assertEquals("Position 7:", "[[ Case3 annotation11 6-9 17-20 28-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 8);
        assertEquals("Position 8:", "[[ Case3 annotation11 6-9 17-20 28-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 9);
        assertEquals("Position 9:", "[[ Case3 annotation11 6-9 17-20 28-29  ]]", currentAnnotations.toString());

         currentAnnotations= getCurrentAnnotations(refIndex, 10);
         assertEquals("Position 10:", "[]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 11);
         assertEquals("Position 11:", "[]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 12);
         assertEquals("Position 12:", "[]", currentAnnotations.toString());

         currentAnnotations= getCurrentAnnotations(refIndex, 13);
        assertEquals("Position 13:", "[[ Case3 annotation8 2-5 13-17 26-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 14);
        assertEquals("Position 14", "[[ Case3 annotation8 2-5 13-17 26-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 15);
        assertEquals("Position 15:", "[[ Case3 annotation8 2-5 13-17 26-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 16);
        assertEquals("Position 16:", "[[ Case3 annotation8 2-5 13-17 26-29  ]]", currentAnnotations.toString());

         currentAnnotations= getCurrentAnnotations(refIndex, 17);
        assertEquals("Position 17:", "[[ Case3 annotation8 2-5 13-17 26-29  ], [ Case3 annotation11 6-9 17-20 28-29  ]]", currentAnnotations.toString());
        
         currentAnnotations= getCurrentAnnotations(refIndex, 18);
        assertEquals("Position 18:", "[[ Case3 annotation11 6-9 17-20 28-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 19);
        assertEquals("Position 19:", "[[ Case3 annotation11 6-9 17-20 28-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 20);
        assertEquals("Position 20:", "[[ Case3 annotation11 6-9 17-20 28-29  ]]", currentAnnotations.toString());

         currentAnnotations= getCurrentAnnotations(refIndex, 21);
        assertEquals("Position 21:", "[]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 22);
        assertEquals("Position 22:", "[]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 23);
        assertEquals("Position 23:", "[]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 24);
        assertEquals("Position 24:", "[]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 25);
        assertEquals("Position 25:", "[]", currentAnnotations.toString());

         currentAnnotations= getCurrentAnnotations(refIndex, 26);
        assertEquals("Position 26:", "[[ Case3 annotation8 2-5 13-17 26-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 27);
        assertEquals("Position 27:", "[[ Case3 annotation8 2-5 13-17 26-29  ]]", currentAnnotations.toString());

         currentAnnotations= getCurrentAnnotations(refIndex, 28);
        assertEquals("Position 28:", "[[ Case3 annotation8 2-5 13-17 26-29  ], [ Case3 annotation11 6-9 17-20 28-29  ]]", currentAnnotations.toString());
         currentAnnotations= getCurrentAnnotations(refIndex, 29);
        assertEquals("Position 29:", "[[ Case3 annotation8 2-5 13-17 26-29  ], [ Case3 annotation11 6-9 17-20 28-29  ]]", currentAnnotations.toString());

    }

    @Test
    public void TestSortedAnnotationsCase4() {

            int refIndex= genome.getReferenceIndex("Chr1");
            ObjectArrayList<Annotation> currentAnnotations;

        currentAnnotations= getCurrentAnnotations(refIndex, 2);
        assertEquals("Position 2:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 3);
        assertEquals("Position 3:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 4);
        assertEquals("Position 4:", "[[ Chr1 annotation12 4-7  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 5);
        assertEquals("Position 5:", "[[ Chr1 annotation12 4-7  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 8);
        assertEquals("Position 8:", "[]", currentAnnotations.toString());

        refIndex= genome.getReferenceIndex("Chr2");

        currentAnnotations= getCurrentAnnotations(refIndex, 2);
        assertEquals("Position 2:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 3);
        assertEquals("Position 3:", "[]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 4);
        assertEquals("Position 4:", "[[ Chr2 annotation13 4-7  ]]", currentAnnotations.toString());

        currentAnnotations= getCurrentAnnotations(refIndex, 5);
        assertEquals("Position 5:", "[[ Chr2 annotation13 4-7  ]]", currentAnnotations.toString());


        currentAnnotations= getCurrentAnnotations(refIndex, 8);
        assertEquals("Position 8:", "[]", currentAnnotations.toString());



    }
}
