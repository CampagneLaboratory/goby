/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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

package samples;import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import it.unimi.dsi.lang.MutableString;

import java.io.FileInputStream;
import java.io.FileNotFoundException;

/**
 * Code sample to demonstrate how to read a compact reads file.
 *
 * @author Fabien Campagne
 *         Date: Jun 24, 2010
 *         Time: 1:15:28 PM
 */
public class ParseReadsCodeSample {
    static {
        try {
            final String inputFilename = "input.compact-reads";
            final MutableString sequence = new MutableString();
            final ReadsReader reader =  new ReadsReader(new FileInputStream(inputFilename));

            for (final Reads.ReadEntry readEntry : reader) {

                ReadsReader.decodeSequence(readEntry, sequence);
                System.out.printf("read-index: %d read-id: %s sequence: %s %n",
                        readEntry.getReadIndex(),
                        readEntry.hasReadIdentifier() ? readEntry.getReadIdentifier() : "",
                        sequence);
            }
        } catch (FileNotFoundException e) {

        }
    }
}
