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
        assertEquals(30, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(50, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(0, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
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
        assertEquals(30, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(50, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(1, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertTrue(samHelper.isReverseStrand());
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
        assertEquals(30, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(3, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(2, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        for (SamSequenceVariation var : vars) {
            System.out.println(var.toString());
        }
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 57, "TTC", 27, "CGT", new byte[]{1, 2, 3}));
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
        assertEquals(30, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(3, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(3, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertTrue(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 52, "TCC", 29, "CTA", new byte[]{3, 2, 1}));
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
        assertEquals(30, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(47, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(3, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(4, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 56, "---", 27, "ATC", new byte[]{1, 2, 3}));
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
        assertEquals(33, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(47, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(3, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(5, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertTrue(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 54, "---", 29, "GAT", new byte[]{3, 2, 1}));
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
        assertEquals(30, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(47, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(3, samHelper.getNumDeletions());
        assertEquals(6, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
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
        assertEquals(30, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(47, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(47, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(3, samHelper.getNumDeletions());
        assertEquals(7, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertTrue(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 52, "TCC", 26, "---", null));
    }

    @Test
    public void testExactStartOfRef() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CTGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        samHelper.setSource(8, sourceRead, sourceQual, "50M", "50", 1, false);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(sourceRead, samHelper.getQuery().toString());
        assertEquals(sourceRead, samHelper.getRef().toString());
        assertEquals(sourceQual, samHelper.getQual().toString());
        assertEquals(0, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(50, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(8, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(0, vars.size());
    }

    @Test
    public void testExactStartOfRefReverse() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CTGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String sourceQual = "XWVUTSRQPONMLKJIHGFEDCBAZYXWVUTSRQPONMLKJIHGFEDCBA";
        samHelper.setSource(9, sourceRead, sourceQual, "50M", "50", 1, true);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(sourceRead, samHelper.getQuery().toString());
        assertEquals(sourceRead, samHelper.getRef().toString());
        assertEquals(sourceQual, samHelper.getQual().toString());
        assertEquals(0, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(50, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(9, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertTrue(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(0, vars.size());
    }

    @Test
    public void testExactEndOfRef() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "GAGTCCTGGCCCCACCCTGTGCTTCCCCTCCGCCTGCTGCACAGGCGCCC";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        samHelper.setSource(10, sourceRead, sourceQual, "50M", "50", 61, false);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(sourceRead, samHelper.getQuery().toString());
        assertEquals(sourceRead, samHelper.getRef().toString());
        assertEquals(sourceQual, samHelper.getQual().toString());
        assertEquals(60, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(50, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(10, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(0, vars.size());
    }

    @Test
    public void testExactEndOfRefReverse() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "GAGTCCTGGCCCCACCCTGTGCTTCCCCTCCGCCTGCTGCACAGGCGCCC";
        final String sourceQual = "XWVUTSRQPONMLKJIHGFEDCBAZYXWVUTSRQPONMLKJIHGFEDCBA";
        samHelper.setSource(11, sourceRead, sourceQual, "50M", "50", 61, true);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(sourceRead, samHelper.getQuery().toString());
        assertEquals(sourceRead, samHelper.getRef().toString());
        assertEquals(sourceQual, samHelper.getQual().toString());
        assertEquals(60, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(50, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(0, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(11, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertTrue(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(0, vars.size());
    }

    @Test
    public void testTrimOneMismatchOne() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "ATGAGGCAGGGCTGGACCCGGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        final String expRead =    "TGAGGCAGGGCTGGACCCGGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String expRef  =    "TGAGGCAGGGCTGGACCCaGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String expQual =    "BCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        samHelper.setSource(12, sourceRead, sourceQual, "1S49M", "18A30", 2, false);
        assertEquals(1, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(expRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(expQual, samHelper.getQual().toString());
        assertEquals(1, samHelper.getPosition());
        assertEquals(1, samHelper.getQueryPosition());
        assertEquals(48, samHelper.getScore());
        assertEquals(49, samHelper.getAlignedLength());
        assertEquals(49, samHelper.getQueryAlignedLength());
        assertEquals(49, samHelper.getTargetAlignedLength());
        assertEquals(1, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(12, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 20, "A", 20, "G", new byte[]{20}));
    }

    @Test
    public void testTrimTwoMismatchOne() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "AAGAGGCAGGGCTGGACCCGGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        final String expRead =    "GAGGCAGGGCTGGACCCGGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String expRef  =    "GAGGCAGGGCTGGACCCaGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String expQual =    "CDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        samHelper.setSource(13, sourceRead, sourceQual, "2S48M", "17A30", 3, false);
        assertEquals(2, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(expRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(expQual, samHelper.getQual().toString());
        assertEquals(2, samHelper.getPosition());
        assertEquals(2, samHelper.getQueryPosition());
        assertEquals(47, samHelper.getScore());
        assertEquals(48, samHelper.getAlignedLength());
        assertEquals(48, samHelper.getQueryAlignedLength());
        assertEquals(48, samHelper.getTargetAlignedLength());
        assertEquals(1, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(13, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 20, "A", 20, "G", new byte[]{20}));
    }

    @Test
    public void testTrimFiveMismatchOne() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "AAAAAGCAGGGCTGGACCCGGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        final String expRead =    "GCAGGGCTGGACCCGGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String expRef  =    "GCAGGGCTGGACCCaGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String expQual =    "FGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        samHelper.setSource(14, sourceRead, sourceQual, "5S45M", "14A30", 6, false);
        assertEquals(5, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(expRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(expQual, samHelper.getQual().toString());
        assertEquals(5, samHelper.getPosition());
        assertEquals(5, samHelper.getQueryPosition());
        assertEquals(44, samHelper.getScore());
        assertEquals(45, samHelper.getAlignedLength());
        assertEquals(45, samHelper.getQueryAlignedLength());
        assertEquals(45, samHelper.getTargetAlignedLength());
        assertEquals(1, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(14, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        for (SamSequenceVariation var : vars) {
            System.out.println(var.toString());
        }
        assertTrue(SamSequenceVariation.contains(vars, 20, "A", 20, "G", new byte[]{20}));
    }

    @Test
    public void testMismatchAtPosTwo() throws IOException {
        final SamHelper samHelper = new SamHelper();
        samHelper.setQualShift((char) -64);
        final String sourceRead = "CAGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String sourceQual = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        final String expRead =    "CAGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String expRef  =    "CtGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC";
        final String expQual =    "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX";
        samHelper.setSource(15, sourceRead, sourceQual, "50M", "1T48", 1, false);
        assertEquals(0, samHelper.getNumLeftClipped());
        assertEquals(0, samHelper.getNumRightClipped());
        assertEquals(expRead, samHelper.getQuery().toString());
        assertEquals(expRef, samHelper.getRef().toString());
        assertEquals(expQual, samHelper.getQual().toString());
        assertEquals(0, samHelper.getPosition());
        assertEquals(0, samHelper.getQueryPosition());
        assertEquals(49, samHelper.getScore());
        assertEquals(50, samHelper.getAlignedLength());
        assertEquals(50, samHelper.getQueryAlignedLength());
        assertEquals(50, samHelper.getTargetAlignedLength());
        assertEquals(1, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(15, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(1, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 2, "T", 2, "A", new byte[]{2}));
    }

/*
    samHelper.setSource(16, "CAAAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC", "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX", "3S47M", "47", 4, false);
    samHelper.setSource(17, "CAAAAACAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTCCC", "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX", "6S44M", "44", 7, false);
    samHelper.setSource(18, "CTGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTCCA", "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX", "49M1S", "49", 1, false);
    samHelper.setSource(19, "CTGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTCAA", "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX", "48M2S", "48", 1, false);
    samHelper.setSource(20, "CTGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCAAAAA", "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX", "45M5S", "45", 1, false);
    samHelper.setSource(21, "CTGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTCAC", "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX", "48M2S", "48", 1, false);
    samHelper.setSource(22, "CTGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTCCTAAC", "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX", "47M3S", "47", 1, false);
    samHelper.setSource(23, "CTGAGGCAGGGCTGGACCCAGTGCCCGCGGCCGCCCTTGCCCTTAAAAAC", "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX", "44M6S", "44", 1, false);
*/
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
        assertEquals(31, samHelper.getPosition());
        assertEquals(1, samHelper.getQueryPosition());
        assertEquals(43, samHelper.getScore());
        assertEquals(48, samHelper.getAlignedLength());
        assertEquals(45, samHelper.getQueryAlignedLength());
        assertEquals(48, samHelper.getTargetAlignedLength());
        assertEquals(2, samHelper.getNumMisMatches());
        assertEquals(0, samHelper.getNumInsertions());
        assertEquals(3, samHelper.getNumDeletions());
        assertEquals(24, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(2, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 41, "CC", 11, "AA", new byte[]{11, 12}));
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
        assertEquals(31, samHelper.getPosition());
        assertEquals(1, samHelper.getQueryPosition());
        assertEquals(43, samHelper.getScore());
        assertEquals(48, samHelper.getAlignedLength());
        assertEquals(48, samHelper.getQueryAlignedLength());
        assertEquals(45, samHelper.getTargetAlignedLength());
        assertEquals(2, samHelper.getNumMisMatches());
        assertEquals(3, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(25, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertFalse(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(2, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 41, "CC", 11, "AA", new byte[]{11, 12}));
        assertTrue(SamSequenceVariation.contains(vars, 56, "---", 27, "ATC", new byte[]{1, 2, 3}));
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
        assertEquals(31, samHelper.getPosition());
        assertEquals(1, samHelper.getQueryPosition());
        assertEquals(43, samHelper.getScore());
        assertEquals(48, samHelper.getAlignedLength());
        assertEquals(48, samHelper.getQueryAlignedLength());
        assertEquals(45, samHelper.getTargetAlignedLength());
        assertEquals(2, samHelper.getNumMisMatches());
        assertEquals(3, samHelper.getNumInsertions());
        assertEquals(0, samHelper.getNumDeletions());
        assertEquals(26, samHelper.getQueryIndex());
        assertEquals(sourceRead.length(), samHelper.getQueryLength());
        assertTrue(samHelper.isReverseStrand());
        List<SamSequenceVariation> vars = samHelper.getSequenceVariations();
        assertEquals(2, vars.size());
        assertTrue(SamSequenceVariation.contains(vars, 41, "CC", 40, "AA", new byte[]{14, 13}));
        assertTrue(SamSequenceVariation.contains(vars, 56, "---", 24, "ATC", new byte[]{24, 23, 22}));
    }

}
