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

package samples;

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Code sample to demonstrate how to iterate through a compact alignment entries file.
 *
 * @author Fabien Campagne
 *         Date: Jun 24, 2010
 */
public class ParseAlignmentEntriesWithIndexCodeSample {
    static {
        try {
            final String inputFilename = "input.entries";

            final AlignmentReader reader = new AlignmentReaderImpl(inputFilename);
            reader.readHeader();
            // we obtain the target index corresponding to chromosome X:
            final int targetIndex=reader.getTargetIdentifiers().get("X");
            final int position=2323122;
            // the following line with use the index (if available), to skip to a genomic position identified by
            // target sequence index and position:
            reader.skipTo(targetIndex,position);
            for (final Alignments.AlignmentEntry alignmentEntry : reader) {

                System.out.printf("query-index: %d target-index: %d score: %f %n",
                        alignmentEntry.getQueryIndex(),
                        alignmentEntry.getTargetIndex(),
                        alignmentEntry.getScore());
            }
        } catch (FileNotFoundException e) {

        } catch (IOException e) {

        }
    }
}
