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

import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.util.motifs.MotifMatcher;
import edu.cornell.med.icb.goby.util.motifs.SubSequenceMotifMatcher;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.log4j.Logger;

import java.io.IOException;

/**
 * Quick took to calculate the frequency of each base in a reads file.
 *
 * @author Fabien Campagne
 *         Date: 6/8/11
 *         Time: 5:31 PM
 */
public class BaseStats {
    private static final Logger LOG = Logger.getLogger(BaseStats.class);

    public static void main(String args[]) throws IOException {

        String filename = args[0];
        MutableString sequence = new MutableString();
        ReadsReader reader = new ReadsReader(filename);
// an array of longs that will hold the count in elements 'A' 'C' 'G' T'
        long[] tallies = new long['T' + 1];
        ProgressLogger progress = new ProgressLogger(LOG);
        progress.start("Starting to analyze base stats.");
        progress.itemsName = "reads";
        // stores frequencies of the CpX motifs
        long cpFreq[] = new long[4];

        final MotifMatcher matchers[] = new MotifMatcher[]{
                new SubSequenceMotifMatcher("CG"),
                new SubSequenceMotifMatcher("CA"),
                new SubSequenceMotifMatcher("CT"),
                new SubSequenceMotifMatcher("CC")
        };

        for (Reads.ReadEntry entry : reader) {
            for (final MotifMatcher matcher : matchers) {
                matcher.newSequence();
            }
            ReadsReader.decodeSequence(entry, sequence);
            for (int i = 0; i < sequence.length(); i++) {
                char c = sequence.charAt(i);
                tallies[c]++;
                int index = 0;
                for (final MotifMatcher matcher : matchers) {
                    matcher.accept(c);
                    if (matcher.matched()) {
                        ++cpFreq[index];
                    }
                    index++;
                }
            }

            progress.lightUpdate();

        }
        progress.done();

        long sum = tallies['A'] + tallies['C'] + tallies['T'] + tallies['G'];
        System.out.printf("tally of bases: %n" +
                "A: %d %g%% %n" +
                "C: %d %g%% %n" +
                "T: %d %g%% %n" +
                "G: %d %g%% %n",
                tallies['A'], percent(tallies['A'], sum),
                tallies['C'], percent(tallies['C'], sum),
                tallies['T'], percent(tallies['T'], sum),
                tallies['G'], percent(tallies['G'], sum)
        );
        long sumcpX = cpFreq[0] + cpFreq[1] + cpFreq[2] + cpFreq[3];
        System.out.printf("tally of CpX motifs: %n" +

                "CpG: %d %g%% %n" +
                "CpA: %d %g%% %n" +
                "CpT: %d %g%% %n" +
                "CpC: %d %g%% %n" +
                "CpT+CpC: %d %g%% %n",
                cpFreq[0], percent(cpFreq[0], sumcpX),
                cpFreq[1], percent(cpFreq[1], sumcpX),
                cpFreq[2], percent(cpFreq[2], sumcpX),
                cpFreq[3], percent(cpFreq[3], sumcpX),
                cpFreq[2]+cpFreq[3], percent(cpFreq[2]+cpFreq[3], sumcpX)
        );
    }

    private static double percent(double tally, double sum) {
        return tally / sum * 100.0;
    }

}
