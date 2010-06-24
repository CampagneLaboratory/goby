/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import it.unimi.dsi.lang.MutableString;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Code sample to demonstrate how to iterate through a compact alignment entries file.
 *
 * @author Fabien Campagne
 *         Date: Jun 24, 2010
 */
public class ParseAlignmentEntriesCodeSample { 
    static {
        try {
            String inputFilename = "input.entries";

            AlignmentReader reader = new AlignmentReader(inputFilename);

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