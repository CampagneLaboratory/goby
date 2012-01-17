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
import org.junit.*;

import java.io.IOException;
import java.util.concurrent.atomic.AtomicReference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * @author Nyasha Chambwe
 * Date: 1/10/12
 * Time: 2:31 PM
 */

public class TestSortedAnnotations {

    SortedAnnotations testObj;
    Annotation [] annotations;

    String[] sequences = {
                "ATTTACCGG", //CpG
                "TAGATACGGAT",//CpG
                "ACTCTAGACTA", //CpT
                "CATTTTGCAAC", //CpA
                "ATCTATGCCTA", //CpG - negative

        };

    RandomAccessSequenceTestSupport genome =new RandomAccessSequenceTestSupport(sequences);

    @Before
    public void setUp() throws IOException {
      /* 1       -1      ENSG00000209346 ENSE00001501617 556070  556142
        1       -1      ENSG00000209349 ENSE00001501620 556239  556304
        1       -1      ENSG00000209350 ENSE00001501621 557859  557930
        2       1       ENSG00000079263 ENSE00000843846 230843545       230843598
        8       -1      ENSG00000077782 ENSE00001289684 38390593        38390698
        8       -1      ENSG00000077782 ENSE00001527273 38390827        38390964*/


      //1       1       ENSG00000209342 ENSE00001501613 554815  554882
      final Annotation annot1 = new Annotation("ENSG00000209342","no-name","1");
      annot1.addSegment(new Segment(554815, 554882, "ENSE00001501613", "+"));

      // 1       1       ENSG00000217866 ENSE00001557252 554882  555926
      final Annotation annot2 = new Annotation("ENSG00000217866", "no-name", "1");
      annot2.addSegment(new Segment(554882, 555926, "ENSE00001557252", "1"));

      //1       1       ENSG00000209343 ENSE00001501614 555725  555900
      Annotation annot3 = new Annotation("ENSG00000209343", "no-name", "1");
      annot3.addSegment(new Segment(555725, 555900, "ENSE00001501614", "1"));

      annotations= new Annotation[]{annot1, annot2, annot3};
    }

    @Test
    public void TestSortedAnnotations(){

        SortedAnnotations testObj= new SortedAnnotations();
        testObj.setAnnotations(annotations);
        testObj.setGenome(genome);


       boolean positiveCase=   testObj.hasOverlappingAnnotations("no-name", 554882);
       assertTrue(positiveCase);
       assertEquals(testObj.currentAnnotations().size(), 2);

       // order: ENSG00000209342, ENSG00000217866
       Annotation first= testObj.currentAnnotations().peek(0);
            
       assertEquals("ENSG00000217866", first.getId());

       Annotation second= testObj.currentAnnotations().peek(1);
       assertEquals("ENSG00000209342", second.getId());
        
       boolean negativeCase=   testObj.hasOverlappingAnnotations("5", 554882);
       assertFalse(negativeCase);

       boolean anotherNegativeCase=   testObj.hasOverlappingAnnotations("no-name", 5882);
       assertFalse(anotherNegativeCase);
        
        
    }




}
