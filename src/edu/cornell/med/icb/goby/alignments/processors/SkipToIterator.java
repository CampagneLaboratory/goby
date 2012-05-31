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

package edu.cornell.med.icb.goby.alignments.processors;

import edu.cornell.med.icb.goby.alignments.Alignments;

import java.io.IOException;

/**
 * Provide a way to add a skipTo interface around a variety of inputs.
 *
 * @author Fabien Campagne
 *         Date: May 1, 2011
 *         Time: 12:12:55 PM
 */
abstract class SkipToIterator {
    /**
     * Return entries from the input reader as if skipTo was called on the iterator input .
     *
     * @param targetIndex Index of the target sequence.
     * @param position    position in the target sequence.
     * @return the next entry that would be returned by skipTo on the reader.
     * @throws IOException If a problem occures obtaining the next entry.
     */
    public abstract Alignments.AlignmentEntry skipTo(int targetIndex, int position) throws IOException;
}
