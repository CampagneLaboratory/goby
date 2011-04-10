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

package edu.cornell.med.icb.goby.alignments;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Apr 9, 2011
 *         Time: 4:29:29 PM
 */
public interface AlignmentReaderFactory {

    /**
     * Returns an implementation of AlignmentReader.
     *
     * @param basename Basename of the alignment.
     * @return implementation of alignment reader, typically customized for some task.
     * @throws java.io.IOException If an error occurs opening a reader.
     */
    public AlignmentReader createReader(String basename) throws IOException;

    /**
     * Returns an array to hold AlignmentReader implementations.
     *
     * @param numElements Number of elements in the array.
     * @return an array of the alignment reader type.
     * @throws java.io.IOException If an error occurs opening a reader.
     */
    AlignmentReader[] createReaderArray(int numElements) throws IOException;

    AlignmentReader createReader(String basename,
                                 int startReferenceIndex, int startPosition,
                                 int endReferenceIndex, int endPosition) throws IOException;
}
