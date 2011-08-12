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

package edu.cornell.med.icb.goby.stats;

import org.junit.Test;
import static org.junit.Assert.assertEquals;

import java.io.StringWriter;
import java.util.Collections;
import java.util.Arrays;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: Apr 7, 2011
 *         Time: 10:14:34 AM
 */
public class TestVCFWriter {
    @Test
    public void testCodeGenotypes() {
        VCFWriter writer = new VCFWriter(new StringWriter());

        assertExpectedGenotype(writer, new String[]{"A", "C"}, 'A', "C", "0/1");
        assertExpectedGenotype(writer, new String[]{"C", "A"}, 'A', "C", "0/1");

        assertExpectedGenotype(writer, new String[]{"C", "A", "T"}, 'A', "C,T", "0/1/2");
        assertExpectedGenotype(writer, new String[]{"C", "A", "T"}, 'C', "A,T", "0/1/2");

    }


    private void assertExpectedGenotype(VCFWriter writer, String[] alleles, char refBase, String altBases, String expectedGenotype) {
        final ObjectArrayList<String> refAlleles = new ObjectArrayList<String>();
        final ObjectArrayList<String> altAlleles = new ObjectArrayList<String>();

        refAlleles.clear();
        refAlleles.add(Character.toString(refBase));
        altAlleles.clear();
        altAlleles.addAll(Arrays.asList(altBases.split(",")));
        final MutableString calculatedGenotype = writer.codeGenotype(alleles, refAlleles, altAlleles);
        assertEquals("Coded genotype must match expected", expectedGenotype, calculatedGenotype.toString());
    }

}
