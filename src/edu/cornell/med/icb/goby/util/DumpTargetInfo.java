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

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: 11/25/11
 *         Time: 3:39 PM
 */
public class DumpTargetInfo {


    public static void main(String args[]) throws IOException {
        String filename = args[0];
        AlignmentReader reader = new AlignmentReaderImpl(filename);
        reader.readHeader();
        IndexedIdentifier targetIds = reader.getTargetIdentifiers();
        DoubleIndexedIdentifier reverse=new DoubleIndexedIdentifier(targetIds);
        for (int targetIndex = 0; targetIndex < targetIds.size(); targetIndex++) {
            System.out.printf("targetIndex=%d id=%s%n", targetIndex, reverse.getId(targetIndex));
        }
        System.out.flush();
        System.exit(0);
    }
}
