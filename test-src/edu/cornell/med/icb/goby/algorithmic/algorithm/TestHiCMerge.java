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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.AlignmentToTextWriter;
import it.unimi.dsi.lang.MutableString;
import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 5/8/12
 *         Time: 4:11 PM
 */
public class TestHiCMerge {

    private static final String EXPECTED_MERGED_INPUT_0 = "{query_index: 537\n" +
            "target_index: 0\n" +
            "position: 165629044\n" +
            "score: 34.0\n" +
            "matching_reverse_strand: true\n" +
            "multiplicity: 1\n" +
            "number_of_mismatches: 2\n" +
            "number_of_indels: 0\n" +
            "query_length: 36\n" +
            "query_aligned_length: 36\n" +
            "target_aligned_length: 36\n" +
            "sequence_variations {\n" +
            "  to: \"C\"\n" +
            "  from: \"T\"\n" +
            "  position: 1\n" +
            "  read_index: 37\n" +
            "}\n" +
            "sequence_variations {\n" +
            "  to: \"G\"\n" +
            "  from: \"T\"\n" +
            "  position: 26\n" +
            "  read_index: 12\n" +
            "}\n" +
            "pair_flags: 161\n" +
            "pair_alignment_link {\n" +
            "  target_index: 0\n" +
            "  position: 165629044\n" +
            "  fragment_index: 1\n" +
            "}\n" +
            "fragment_index: 0\n" +
            "}\n" +
            "{query_index: 537\n" +
            "target_index: 0\n" +
            "position: 165629044\n" +
            "score: 34.0\n" +
            "matching_reverse_strand: true\n" +
            "multiplicity: 1\n" +
            "number_of_mismatches: 2\n" +
            "number_of_indels: 0\n" +
            "query_length: 36\n" +
            "query_aligned_length: 36\n" +
            "target_aligned_length: 36\n" +
            "sequence_variations {\n" +
            "  to: \"C\"\n" +
            "  from: \"T\"\n" +
            "  position: 1\n" +
            "  read_index: 37\n" +
            "}\n" +
            "sequence_variations {\n" +
            "  to: \"G\"\n" +
            "  from: \"T\"\n" +
            "  position: 26\n" +
            "  read_index: 12\n" +
            "}\n" +
            "pair_alignment_link {\n" +
            "  target_index: 0\n" +
            "  position: 165629044\n" +
            "  fragment_index: 0\n" +
            "}\n" +
            "fragment_index: 1\n" +
            "}\n" +
            "Set smallestSplitQueryIndex=0\n" +
            "Set largestSplitQueryIndex=537\n" +
            "Set targetIdentifiers: {0->1\n" +
            "}\n" +
            "Set targetLengths: {0->249250621\n" +
            "}\n" +
            "Set numQueries=22819628\n" +
            "Set numTargets=1\n" +
            "Set statistics { {part1.number.aligned.reads=1, part1.basename.full=mantis-1355, part2.basename.full=mantis-1355, part2.basename=mantis-1355, part1.max.query.index=537, part1.number.of.queries=22819628, part2.max.query.index=537, part2.number-of-entries-written=1, part2.number.aligned.reads=1, part1.number-of-entries-written=1, part2.min.query.index=0, part1.min.query.index=0, part2.number.of.queries=22819628, part1.basename=mantis-1355} }\n" +
            "Closed\n";


    @Test
    public void mergeHiCTestCase1() throws IOException {

        AlignmentToTextWriter writer = new AlignmentToTextWriter();
        HiCMerge merger = new HiCMerge();

        final String basenameA = "test-data/alignments/mantis-1355/mantis-1355.entries";
        AlignmentReader readerA = new AlignmentReaderImpl(basenameA);
        AlignmentReader readerB = new AlignmentReaderImpl(basenameA);
        merger.merge(readerA, readerB, writer);
        merger.transferHeader(basenameA,basenameA,writer);
        writer.close();

       // System.out.println("--");
       // System.out.println(writer.getTextOutput());
       // System.out.println("--");
        assertEquals( "The text must match",EXPECTED_MERGED_INPUT_0, writer.getTextOutput().toString());


    }


}
