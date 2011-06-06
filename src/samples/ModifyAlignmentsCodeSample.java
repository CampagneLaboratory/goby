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

package samples;

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.IterateAlignments;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Jun 24, 2010
 *         Time: 2:43:44 PM
 */
public class ModifyAlignmentsCodeSample {
    static {

        try {
            final IterateAlignments iterator = new IterateAlignments() {

                public void processAlignmentEntry(final AlignmentReader alignmentReader, final Alignments.AlignmentEntry alignmentEntry) {
                  /*
                    int tagValue =???

                    final Alignments.AlignmentEntry revisedEntry = Alignments.AlignmentEntry.newBuilder()
                            .mergeFrom(alignmentEntry).setMyTagValue(tagValue).build();

                    // do something with  revisedEntry
                    */
                }
            };
            final String[] basenames = {"input1.entries", "input2.entries"};

            iterator.parseIncludeReferenceArgument("1,2,X");
            iterator.iterate(basenames);


        } catch (IOException e) {

        }
    }
}
