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

package edu.cornell.med.icb.goby.util;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;

/**
 * @author Fabien Campagne
 *         Date: Apr 19, 2011
 *         Time: 9:22:22 PM
 */
public class TestSimulateBisulfiteReads {

    @Test
    public void testBisulfiteForwardStrand() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{1, 0, 1, 0, 1, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(true);
        simulator.setDoReverseStrand(false);
        simulator.setBisulfiteTreatment(true);
        simulator.process(rates, "CCCCCCCCC", 0, stringBuffer);
        String expected=
                "@0 reference: null startPosition: 0 strand: +1 met: 1 read-index: 1 met: 3 read-index: 3 met: 5 read-index: 5 met: 7 read-index: 7 met: 9 read-index: 9  1 2 3 4 5 6 7 8 9 \n" +
                "CTCTCTCTC\n" +
                "+\n" +
                "hhhhhhhhh\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

     @Test
    public void testBisulfiteReverseStrand() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{0,0,1,  0, 0, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(false);
        simulator.setDoReverseStrand(true);
        simulator.setBisulfiteTreatment(true);
        simulator.process(rates, "ACTGGG", 0, stringBuffer);
        String expected=
                "@0 reference: null startPosition: 0 strand: -1 met: 3 read-index: 4  1 2 3 4 5 6 \n" +
                "TTCAGT\n" +
                "+\n" +
                "hhhhhh\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

      @Test
    public void testMutateForwardStrand() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{1,0,  0, 0, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(true);
        simulator.setDoReverseStrand(false);
        simulator.setBisulfiteTreatment(false);
        simulator.process(rates, "ACTCGG", 0, stringBuffer);
        String expected=
                "@0 reference: null startPosition: 0 strand: +1 mut: 2 read-index: 2  1 2 3 4 5 6 \n" +
                "AGTCGG\n" +
                "+\n" +
                "hhhhhh\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

      @Test
    public void testMutateReverseStrand() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{0,1,0,  0, 0, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(false);
        simulator.setDoReverseStrand(true);
        simulator.setBisulfiteTreatment(false);
        simulator.process(rates, "ACTCGG", 0, stringBuffer);
        String expected=
                "@0 reference: null startPosition: 0 strand: -1 mut: 2 read-index: 5  1 2 3 4 5 6 \n" +
                "CGGAGT\n" +
                "+\n" +
                "hhhhhh\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }
}
