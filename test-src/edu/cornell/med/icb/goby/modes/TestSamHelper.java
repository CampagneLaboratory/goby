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

package edu.cornell.med.icb.goby.modes;

import org.junit.Test;
import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

/**
 * Test construction of ref via sam data.
 */
public class TestSamHelper {

    @Test
    public void testSimple() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CCGCCCTTGCCCTTCCTCCCTTCCCTTTCGGAGTCCTGGCCCCACCCTGT";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        samHelper.setSource(0, sourceRead, sourceQual, "50M", "50", 31, false);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(sourceRead, samHelper.getQuery().toString());
        assertEquals(sourceRead, samHelper.getRef().toString());
        assertEquals(sourceQual, samHelper.getQual().toString());
        assertEquals(50, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(0, vars.size());
    }

    @Test
    public void testSimpleReverse() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CCGCCCTTGCCCTTCCTCCCTTCCCTTTCGGAGTCCTGGCCCCACCCTGT";
        final String sourceQual = "XWVUTSRQPONMLKJIHGFEDCBAZYXWVUTSRQPONMLKJIHGFEDCBA";
        samHelper.setSource(1, sourceRead, sourceQual, "50M", "50", 31, true);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(sourceRead, samHelper.getQuery().toString());
        assertEquals(sourceRead, samHelper.getRef().toString());
        assertEquals(sourceQual, samHelper.getQual().toString());
        assertEquals(50, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(0, vars.size());
    }

    @Test
    public void testMismatches() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CCGCCCTTGCCCTTCCTCCCTTCCCTCGTGGAGTCCTGGCCCCACCCTGT";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        final String expRef =     "CCGCCCTTGCCCTTCCTCCCTTCCCTttcGGAGTCCTGGCCCCACCCTGT";
        samHelper.setSource(2, sourceRead, sourceQual, "50M", "26TTC21", 31, false);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(sourceRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(sourceQual, samHelper.getQual().toString());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 57, "TTC", 27, "CGT", new int[]{1, 2, 3}));
    }

    @Test
    public void testMismatchesReverse() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CCGCCCTTGCCCTTCCTCCCTCTACTTTCGGAGTCCTGGCCCCACCCTGT";
        final String sourceQual = "XWVUTSRQPONMLKJIHGFEDCBAZYXWVUTSRQPONMLKJIHGFEDCBA";
        final String expRef =     "CCGCCCTTGCCCTTCCTCCCTtccCTTTCGGAGTCCTGGCCCCACCCTGT";
        samHelper.setSource(3, sourceRead, sourceQual, "50M", "21TCC26", 31, true);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(sourceRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(sourceQual, samHelper.getQual().toString());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 52, "TCC", 29, "CTA", new int[]{3, 2, 1}));
    }

    @Test
    public void testInsert() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CCGCCCTTGCCCTTCCTCCCTTCCCTATCTTCGGAGTCCTGGCCCCACCC";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        final String expRef =     "CCGCCCTTGCCCTTCCTCCCTTCCCT---TTCGGAGTCCTGGCCCCACCC";
        samHelper.setSource(4, sourceRead, sourceQual, "26M3I21M", "47", 31, false);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(sourceRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(sourceQual, samHelper.getQual().toString());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(47, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 56, "---", 27, "ATC", new int[]{1, 2, 3}));
    }

    @Test
    public void testInsertReverse() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CCCTTGCCCTTCCTCCCTTCCGATCTTTCGGAGTCCTGGCCCCACCCTGT";
        final String sourceQual = "XWVUTSRQPONMLKJIHGFEDCBAZYXWVUTSRQPONMLKJIHGFEDCBA";
        final String expRef =     "CCCTTGCCCTTCCTCCCTTCC---CTTTCGGAGTCCTGGCCCCACCCTGT";
        samHelper.setSource(5, sourceRead, sourceQual, "21M3I26M", "47", 34, true);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(sourceRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(sourceQual, samHelper.getQual().toString());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(47, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 54, "---", 29, "GAT", new int[]{3, 2, 1}));
    }

    @Test
    public void testDelete() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CCGCCCTTGCCCTTCCTCCCTTCCCTGGAGTCCTGGCCCCACCCTGT";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTU";
        final String expRef =     "CCGCCCTTGCCCTTCCTCCCTTCCCTTTCGGAGTCCTGGCCCCACCCTGT";
        final String expRead =    "CCGCCCTTGCCCTTCCTCCCTTCCCT---GGAGTCCTGGCCCCACCCTGT";
        final String expQual =    "ABCDEFGHIJKLMNOPQRSTUVWXYZ___ABCDEFGHIJKLMNOPQRSTU";
        samHelper.setMinQualValue('_');
        samHelper.setSource(6, sourceRead, sourceQual, "26M3D21M", "26^TTC21", 31, false);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(expRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(expQual, samHelper.getQual().toString());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(47, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 57, "TTC", 26, "---", null));
    }

    @Test
    public void testDeleteReverse() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CCGCCCTTGCCCTTCCTCCCTCTTTCGGAGTCCTGGCCCCACCCTGT";
        final String sourceQual = "UTSRQPONMLKJIHGFEDCBAZYXWVUTSRQPONMLKJIHGFEDCBA";
        final String expRef =     "CCGCCCTTGCCCTTCCTCCCTTCCCTTTCGGAGTCCTGGCCCCACCCTGT";
        final String expRead =    "CCGCCCTTGCCCTTCCTCCCT---CTTTCGGAGTCCTGGCCCCACCCTGT";
        final String expQual =    "UTSRQPONMLKJIHGFEDCBA___ZYXWVUTSRQPONMLKJIHGFEDCBA";
        samHelper.setMinQualValue('_');
        samHelper.setSource(7, sourceRead, sourceQual, "21M3D26M", "21^TCC26", 31, true);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(expRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(expQual, samHelper.getQual().toString());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(47, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 52, "TCC", 26, "---", null));
    }

    @Test
    public void testTrimMismatchAndDelete() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "ACGCCCTTGCAATTCCTCCCTTCCCTTTCGGCCTGGCCCCACCCTGA";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTU";
        final String expRead =    "CGCCCTTGCAATTCCTCCCTTCCCTTTCGG---CCTGGCCCCACCCTG";
        final String expRef  =    "CGCCCTTGCccTTCCTCCCTTCCCTTTCGGAGTCCTGGCCCCACCCTG";
        final String expQual =    "BCDEFGHIJKLMNOPQRSTUVWXYZABCDE___FGHIJKLMNOPQRST";
        samHelper.setMinQualValue('_');
        samHelper.setSource(24, sourceRead, sourceQual, "1S30M3D15M1S", "9CC19^AGT15", 32, false);
        assertEquals(1, samHelper.getNumLeftClipped());
        assertEquals(1, samHelper.getNumRightClipped());
        assertEquals(expRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(expQual, samHelper.getQual().toString());
        assertEquals(43, samHelper.getScore());
        assertEquals(48, samHelper.getAlignedLength());
        assertEquals(45, samHelper.getQueryAlignedLength());
        assertEquals(48, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(2, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 41, "CC", 11, "AA", new int[]{11, 12}));
        assertTrue(SamSequenceVariation.contains(vars, 62, "AGT", 31, "---", null));
    }

    @Test
    public void testTrimMismatchAndInsert() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "ACGCCCTTGCAATTCCTCCCTTCCCTATCTTCGGAGTCCTGGCCCCACCA";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        final String expRead =    "CGCCCTTGCAATTCCTCCCTTCCCTATCTTCGGAGTCCTGGCCCCACC";
        final String expRef  =    "CGCCCTTGCccTTCCTCCCTTCCCT---TTCGGAGTCCTGGCCCCACC";
        final String expQual =    "BCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVW";
        samHelper.setSource(25, sourceRead, sourceQual, "1S25M3I20M1S", "9CC34", 32, false);
        assertEquals(1, samHelper.getNumLeftClipped());
        assertEquals(1, samHelper.getNumRightClipped());
        assertEquals(expRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(expQual, samHelper.getQual().toString());
        assertEquals(43, samHelper.getScore());
        assertEquals(48, samHelper.getAlignedLength());
        assertEquals(48, samHelper.getQueryAlignedLength());
        assertEquals(45, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(2, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 41, "CC", 11, "AA", new int[]{11, 12}));
        assertTrue(SamSequenceVariation.contains(vars, 56, "---", 27, "ATC", new int[]{1, 2, 3}));
    }

    @Test
    public void testTrimMismatchAndInsertReverse() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "ACGCCCTTGCAATTCCTCCCTTCCCTATCTTCGGAGTCCTGGCCCCACCA";
        final String sourceQual = "XWVUTSRQPONMLKJIHGFEDCBAZYXWVUTSRQPONMLKJIHGFEDCBA";
        final String expRead =    "CGCCCTTGCAATTCCTCCCTTCCCTATCTTCGGAGTCCTGGCCCCACC";
        final String expRef  =    "CGCCCTTGCccTTCCTCCCTTCCCT---TTCGGAGTCCTGGCCCCACC";
        final String expQual =    "WVUTSRQPONMLKJIHGFEDCBAZYXWVUTSRQPONMLKJIHGFEDCB";
        samHelper.setSource(26, sourceRead, sourceQual, "1S25M3I20M1S", "9CC34", 32, true);
        assertEquals(1, samHelper.getNumLeftClipped());
        assertEquals(1, samHelper.getNumRightClipped());
        assertEquals(expRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(expQual, samHelper.getQual().toString());
        assertEquals(43, samHelper.getScore());
        assertEquals(48, samHelper.getAlignedLength());
        assertEquals(48, samHelper.getQueryAlignedLength());
        assertEquals(45, samHelper.getTargetAlignedLength());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(2, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 41, "CC", 40, "AA", new int[]{14, 13}));
        assertTrue(SamSequenceVariation.contains(vars, 56, "---", 24, "ATC", new int[]{24, 23, 22}));
    }

}
