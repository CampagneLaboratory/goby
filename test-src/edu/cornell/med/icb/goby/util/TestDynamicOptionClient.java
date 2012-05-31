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

package edu.cornell.med.icb.goby.util;

import edu.cornell.med.icb.goby.stats.AnnotationAveragingWriter;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import org.junit.Test;

import static junit.framework.Assert.assertTrue;

/**
 * @author Fabien Campagne
 *         Date: 2/1/12
 *         Time: 8:32 AM
 */
public class TestDynamicOptionClient {
    private static class MyOptionClass {
        public final static DynamicOptionClient doc = new DynamicOptionClient(MyOptionClass.class, "option1:Some text:", "option2:an int:-3");
    }

    @Test
    public void test1() {
        MyOptionClass testInstance = new MyOptionClass();
        assertTrue("option must be accepted", testInstance.doc.acceptsOption("MyOptionClass:option1:valueOption1"));
        assertTrue("option must be accepted", testInstance.doc.acceptsOption("MyOptionClass:option2=50"));

    }

    @Test
    public void test2(){

        final DynamicOptionClient doc = new DynamicOptionClient(AnnotationAveragingWriter.class,
                "annotations:annotation filename:",
                "write-counts:boolean, when true write C and Cm for regions:false",
                "estimate-intra-group-differences: boolean, true indicates that pair-wise differences for sample in the same group should be tallied and written to the output. False indicates regular output.:false",
                "estimate-empirical-P: boolean, true: activates estimation of the empirical p-value.:false",
                "combinator: string, the method to combine p-values, one of qfast, average, sum, max.:sum",
                "serialized-estimator-filename: string, the path to a serialized version of the density estimator populated with the empirical null-distribution.:"
        );
        assertTrue(doc.acceptsOption("AnnotationAveragingWriter:write-counts=false"));
        assertTrue(doc.acceptsOption("AnnotationAveragingWriter:annotations=filename"));
    }
}
