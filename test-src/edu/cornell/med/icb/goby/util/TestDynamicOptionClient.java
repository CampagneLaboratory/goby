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

import org.junit.Test;

import static junit.framework.Assert.assertTrue;

/**
 * @author Fabien Campagne
 *         Date: 2/1/12
 *         Time: 8:32 AM
 */
public class TestDynamicOptionClient {
    private static class MyOptionClass  {
        public final static DynamicOptionClient doc = new DynamicOptionClient(MyOptionClass.class, "option1:Some text:", "option2:an int:-3");
    }

    @Test
    public void test1() {
        MyOptionClass testInstance = new MyOptionClass();
        assertTrue("option must be accepted",testInstance.doc.acceptsOption("MyOptionClass:option1:valueOption1"));
        assertTrue("option must be accepted",testInstance.doc.acceptsOption("MyOptionClass:option2=50"));

    }
}
