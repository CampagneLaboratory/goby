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
    public void testBisulfiteForwardTrueRates() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{1, 0, 1, 0, 1, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(true);
        simulator.setDoReverseStrand(false);
        simulator.setBisulfiteTreatment(true);
        simulator.writeTrueRatesBisulfite(rates, "CCCCCCCCC", 0, stringBuffer);
        String expected =
                "index\tmethylationRate\tstrand\tchromosome\tposition\tfromBase\ttoBase\tcontext\n" +
                        "0\t1.00000\t+1\tnull\t1\tC\tC\t>C<CCCCCCCC\n" +
                        "1\t0.00000\t+1\tnull\t2\tC\tT\tC>C<CCCCCCC\n" +
                        "2\t1.00000\t+1\tnull\t3\tC\tC\tCC>C<CCCCCC\n" +
                        "3\t0.00000\t+1\tnull\t4\tC\tT\tCCC>C<CCCCC\n" +
                        "4\t1.00000\t+1\tnull\t5\tC\tC\tCCCC>C<CCCC\n" +
                        "5\t0.00000\t+1\tnull\t6\tC\tT\tCCCCC>C<CCC\n" +
                        "6\t1.00000\t+1\tnull\t7\tC\tC\tCCCCCC>C<CC\n" +
                        "7\t0.00000\t+1\tnull\t8\tC\tT\tCCCCCCC>C<C\n" +
                        "8\t1.00000\t+1\tnull\t9\tC\tC\tCCCCCCCC>C<\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

    @Test
    public void testBisulfiteForwardTrueRates2() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{1, 0, 1, 0, 1, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(true);
        simulator.setDoReverseStrand(false);
        simulator.setBisulfiteTreatment(true);
        simulator.writeTrueRatesBisulfite(rates, "ACTGCACT", 0, stringBuffer);
        String expected =
                "index\tmethylationRate\tstrand\tchromosome\tposition\tfromBase\ttoBase\tcontext\n" +
                        "1\t1.00000\t+1\tnull\t2\tC\tC\tA>C<TGCACT\n" +
                        "4\t0.00000\t+1\tnull\t5\tC\tT\tACTG>C<ACT\n" +
                        "6\t1.00000\t+1\tnull\t7\tC\tC\tACTGCA>C<T\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

    @Test
    public void testBisulfiteReverseTrueRates() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{1, 0, 1, 0, 1, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(false);
        simulator.setDoReverseStrand(true);
        simulator.setBisulfiteTreatment(true);
        simulator.writeTrueRatesBisulfite(rates, "GGGGGGGG", 0, stringBuffer);
        String expected =
                "index\tmethylationRate\tstrand\tchromosome\tposition\tfromBase\ttoBase\tcontext\n" +
                        "0\t1.00000\t-1\tnull\t8\tC\tC\tGGGGGGG>G<\n" +
                        "1\t0.00000\t-1\tnull\t7\tC\tT\tGGGGGG>G<G\n" +
                        "2\t1.00000\t-1\tnull\t6\tC\tC\tGGGGG>G<GG\n" +
                        "3\t0.00000\t-1\tnull\t5\tC\tT\tGGGG>G<GGG\n" +
                        "4\t1.00000\t-1\tnull\t4\tC\tC\tGGG>G<GGGG\n" +
                        "5\t0.00000\t-1\tnull\t3\tC\tT\tGG>G<GGGGG\n" +
                        "6\t1.00000\t-1\tnull\t2\tC\tC\tG>G<GGGGGG\n" +
                        "7\t0.00000\t-1\tnull\t1\tC\tT\t>G<GGGGGGG\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

    @Test
    public void testBisulfiteReverseTrueRatesNothing() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{1, 0, 1, 0, 1, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(false);
        simulator.setDoReverseStrand(true);
        simulator.setBisulfiteTreatment(true);
        simulator.writeTrueRatesBisulfite(rates, "CCCCCCCCC", 0, stringBuffer);
        // no G on reverse strand:
        String expected =
                "index\tmethylationRate\tstrand\tchromosome\tposition\tfromBase\ttoBase\tcontext\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

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
        String expected =
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
        DoubleList rates = DoubleArrayList.wrap(new double[]{0, 0, 1, 0, 0, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(false);
        simulator.setDoReverseStrand(true);
        simulator.setBisulfiteTreatment(true);
        // ref: ACTGGG
        // rev: CCCAGT
        // bis: TTCAGT read-index=3
        // rev: ACTGAA pos=4
        // pos: 123456
        simulator.process(rates, "ACTGGG", 0, stringBuffer);
        String expected =
                "@0 reference: null startPosition: 0 strand: -1 met: 4 read-index: 3  1 2 3 4 5 6 \n" +
                        "TTCAGT\n" +
                        "+\n" +
                        "hhhhhh\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

    @Test
    public void testBisulfiteReverseTrueRatesLonger() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{0, 1, 0, 0, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(false);
        simulator.setDoReverseStrand(true);
        simulator.setBisulfiteTreatment(true);
        simulator.writeTrueRatesBisulfite(rates, "TTCGTGTGTGTTTAAAGCTT", 0, stringBuffer);
        String expected =
                "index\tmethylationRate\tstrand\tchromosome\tposition\tfromBase\ttoBase\tcontext\n" +
                        "3\t0.00000\t-1\tnull\t17\tC\tT\tGTGTTTAAA>G<CTT\n" +
                        "10\t1.00000\t-1\tnull\t10\tC\tC\tTTCGTGTGT>G<TTTAAAGCTT\n" +
                        "12\t0.00000\t-1\tnull\t8\tC\tT\tTTCGTGT>G<TGTTTAAAGC\n" +
                        "14\t0.00000\t-1\tnull\t6\tC\tT\tTTCGT>G<TGTGTTTAAA\n" +
                        "16\t0.00000\t-1\tnull\t4\tC\tT\tTTC>G<TGTGTGTTTA\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

    @Test
    public void testBisulfiteReverseStrandFromNotZero() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(5);
        DoubleList rates = DoubleArrayList.wrap(new double[]{0, 1, 0, 0, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(false);
        simulator.setDoReverseStrand(true);
        simulator.setBisulfiteTreatment(true);
        simulator.setReadLength(15);
        // ref: TTCGTGTGT>G<TTTAAAGCTT forward
        // rev: AAGCTTTAAA>C<ACACACGAA ref reverse complemented
        // bis: AAGTTTTAAA>C<ATATATGAA read-index=11
        // rev: TTCATATAT>G<TTTAAAACTT pos=10
        // pos: 123456789 10
        // TTCGTGTGTGTTTAAAGCTT  forward reference
        // AAGCTTTAAACACACACGAA  reversed reference
        simulator.process(rates, "TTCGTGTGTGTTTAAAGCTT", 5, stringBuffer);
        String expected =
                "@0 reference: null startPosition: 0 strand: -1 met: 15 read-index: 6  6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 \n" +
                        "TTAAACATATATGAA\n" +
                        "+\n" +
                        "hhhhhhhhhhhhhhh\n" +
                        "@1 reference: null startPosition: 4 strand: -1 met: 15 read-index: 10  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 \n" +
                        "AGTTTTAAACATATA\n" +
                        "+\n" +
                        "hhhhhhhhhhhhhhh\n" +
                        "@2 reference: null startPosition: 0 strand: -1 met: 15 read-index: 6  6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 \n" +
                        "TTAAACATATATGAA\n" +
                        "+\n" +
                        "hhhhhhhhhhhhhhh\n" +
                        "@3 reference: null startPosition: 2 strand: -1 met: 15 read-index: 8  8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 \n" +
                        "TTTTAAACATATATG\n" +
                        "+\n" +
                        "hhhhhhhhhhhhhhh\n" +
                        "@4 reference: null startPosition: 4 strand: -1 met: 15 read-index: 10  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 \n" +
                        "AGTTTTAAACATATA\n" +
                        "+\n" +
                        "hhhhhhhhhhhhhhh\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

    @Test
    public void testMutateForwardStrand() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{1, 0, 0, 0, 0, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(true);
        simulator.setDoReverseStrand(false);
        simulator.setBisulfiteTreatment(false);
        simulator.process(rates, "ACTCGG", 0, stringBuffer);
        String expected =
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
        DoubleList rates = DoubleArrayList.wrap(new double[]{0, 1, 0, 0, 0, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(false);
        simulator.setDoReverseStrand(true);
        simulator.setBisulfiteTreatment(false);
        //ref: ACTCGG
        //rev: CCGAGT
        //mut: CGGAGT  read-index=2
        //rev: ACTCCG  position=5
        simulator.process(rates, "ACTCGG", 0, stringBuffer);
        String expected =
                "@0 reference: null startPosition: 0 strand: -1 mut: 5 read-index: 2  1 2 3 4 5 6 \n" +
                        "CGGAGT\n" +
                        "+\n" +
                        "hhhhhh\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

@Test
    public void testMutateTrueRates() throws IOException {
        SimulateBisulfiteReads simulator = new SimulateBisulfiteReads();
        simulator.setNumRepeats(1);
        DoubleList rates = DoubleArrayList.wrap(new double[]{1, 0, 1, 0, 1, 0});
        StringWriter stringBuffer = new StringWriter(1000);
        simulator.setDoForwardStrand(true);
        simulator.setDoReverseStrand(false);
        simulator.setBisulfiteTreatment(false);
        simulator.writeTrueRatesMutations(rates, "CCCCCCCCC", 0, stringBuffer);
        String expected =
                "index\tmethylationRate\tstrand\tchromosome\tposition\tfromBase\ttoBase\tcontext\n" +
                        "0\t1.00000\t+1\tnull\t1\tC\tG\t>C<CCCCCCCC\n" +
                        "1\t0.00000\t+1\tnull\t2\tC\tG\tC>C<CCCCCCC\n" +
                        "2\t1.00000\t+1\tnull\t3\tC\tG\tCC>C<CCCCCC\n" +
                        "3\t0.00000\t+1\tnull\t4\tC\tG\tCCC>C<CCCCC\n" +
                        "4\t1.00000\t+1\tnull\t5\tC\tG\tCCCC>C<CCCC\n" +
                        "5\t0.00000\t+1\tnull\t6\tC\tG\tCCCCC>C<CCC\n" +
                        "6\t1.00000\t+1\tnull\t7\tC\tG\tCCCCCC>C<CC\n" +
                        "7\t0.00000\t+1\tnull\t8\tC\tG\tCCCCCCC>C<C\n" +
                        "8\t1.00000\t+1\tnull\t9\tC\tG\tCCCCCCCC>C<\n";
        assertEquals(expected, stringBuffer.getBuffer().toString());
    }

}
