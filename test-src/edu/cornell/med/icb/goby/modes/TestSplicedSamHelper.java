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

package edu.cornell.med.icb.goby.modes;

import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

/**
 * @author Fabien Campagne
 *         Date: 3/29/12
 *         Time: 12:11 PM
 */
public class TestSplicedSamHelper {
    /*
9068319	0	1	18339	2	28M6371N7M	*	0	0	CCTGCACCTGGCTCCGGCTCTGCTCTACCTGCTGA	aaa^`a``a^^aaaa``aaa_V__`_X]`a`aa_[	MD:Z:35	NH:i:3	NM:i:0	SM:i:2	XQ:i:40	X2:i:40	XS:A:-
11090122	0	1	18345	1	22M6371N13M	*	0	0	CCTGGCTCCGGCTCTGCTCTACCTGCTGAAGATGT	Xa^`U\``]]Y`ZZ\ZYZ\\\Z`ZQ\XJO\VGOQQ	MD:Z:35	NH:i:4	NM:i:0	SM:i:1	XQ:i:40	X2:i:40	XS:A:-
12986491	0	1	18345	2	22M6371N13M	*	0	0	CCTGGCTCCGGCTCTGCTCTACCTGCTGCAGATGT	__a`a_^__^`^`^^^^^^]V[``^\YTFV\XXYS	MD:Z:28A6	NH:i:3	NM:i:1	SM:i:2	XQ:i:40	X2:i:40	XS:A:-

     */

    @Test
    public void testLimits() {
        final SplicedSamHelper samHelper = new SplicedSamHelper();

        SplicedSamHelper.Limits[] limits = samHelper.getLimits(18339, "28M6371N7M", "28A6");
        assertEquals(2, limits.length);
        assertEquals(18339, limits[0].position);
        assertEquals(18339 + 28 + 6371, limits[1].position);

        assertEquals(0,  limits[0].cigarStart);
        assertEquals(3,  limits[0].cigarEnd);
        assertEquals(8,  limits[1].cigarStart);
        assertEquals(10, limits[1].cigarEnd);

        assertEquals(0,  limits[0].readStart);
        assertEquals(28, limits[0].readEnd);
        assertEquals(28, limits[1].readStart);
        assertEquals(35, limits[1].readEnd);

        assertEquals("28", limits[0].md);
        assertEquals("A6", limits[1].md);



    }

    @Test
    public void testSpliced() throws IOException {
        final SplicedSamHelper samHelper = new SplicedSamHelper();
        String bases_0_28 = "CCTGCACCTGGCTCCGGCTCTGCTCTAC";
        String quals_0_28 = "aaa^`a``a^^aaaa``aaa_V__`_X]";
        final String bases_28_35 = "CTGCTGA";
        final String sourceRead = bases_0_28 + bases_28_35;
        final String quals_28_35 = "`a`aa_[";
        final String sourceQual = quals_0_28 + quals_28_35;
        samHelper.setSource(0, sourceRead, sourceQual, "28M6371N7M", "35", 18339, false);
        assertEquals(2, samHelper.getNumEntries());
        samHelper.setEntryCursor(0);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(bases_0_28, samHelper.getQuery().toString());
        assertEquals(bases_0_28, samHelper.getRef().toString());
        assertEquals(quals_0_28, samHelper.getQual().toString());
        assertEquals(18339, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(28, samHelper.getScore());
        assertEquals(28, samHelper.getAlignedLength());
        assertEquals(28, samHelper.getQueryAlignedLength());
        assertEquals(28, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(0, samHelper.getQueryIndex());
        assertEquals(bases_0_28.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(0, vars.size());

        samHelper.setEntryCursor(1);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(bases_28_35, samHelper.getQuery().toString());
        assertEquals(bases_28_35, samHelper.getRef().toString());
        assertEquals(quals_28_35, samHelper.getQual().toString());
        assertEquals(18339+6371+bases_0_28.length(), samHelper.getPosition());
    }

}
