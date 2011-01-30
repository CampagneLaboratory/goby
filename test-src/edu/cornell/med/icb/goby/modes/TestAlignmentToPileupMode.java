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

import edu.cornell.med.icb.goby.alignments.AlignmentWriter;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import static org.junit.Assert.fail;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: Jan 29, 2011
 *         Time: 4:57:59 PM
 */
public class TestAlignmentToPileupMode {

    final String basename1 = "test-data/alignments/pileup/1";

    @Test
    public void testMutation() {

        final StringWriter output = new StringWriter();

        AlignmentToPileupMode pileup = new AlignmentToPileupMode();
        pileup.setInputFilename(basename1);
        pileup.setOutputWriter(new PrintWriter(output));
        pileup.setFormat(AlignmentToPileupMode.OutputFormat.ONE_PER_LINE);
        pileup.initializeIterator("1,0","1,50");
 //       pileup.setGenomeFilename("test-data/alignments/pileup/genome.fasta");
        try {

            pileup.execute();
            System.out.println(output.toString());

        } catch (IOException e) {
            fail("An exception was caught" + e.getMessage());
        }
    }

    @Before
    public void setUp() throws IOException {


        AlignmentWriter writer = new AlignmentWriter(basename1);
        final Alignments.AlignmentEntry.Builder builder = Alignments.AlignmentEntry.newBuilder();
        builder.setQueryLength(11);
        builder.setQueryIndex(1);
        builder.setTargetIndex(0);
        builder.setPosition(10);
        builder.setQueryAlignedLength(10);
        final Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder();
        varBuilder.setFrom("A");  // A at pos 15
        varBuilder.setTo("T");
        varBuilder.setPosition(3);
        varBuilder.setReadIndex(9);
        builder.addSequenceVariations(varBuilder.build());
        builder.setMatchingReverseStrand(true);
        final Alignments.AlignmentEntry entry = builder.build();
        writer.appendEntry(entry);
        IndexedIdentifier ids=new IndexedIdentifier();
        ids.put(new MutableString("1"),0);
        writer.setTargetIdentifiers(ids);
        int[] targetLengths = new int[]{50};
        writer.setTargetLengths(targetLengths);
        writer.setSorted(true);
        writer.close();
    }


}
