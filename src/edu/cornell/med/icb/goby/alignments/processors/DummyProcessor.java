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

import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.alignments.processors.AlignmentProcessorInterface;
import edu.cornell.med.icb.goby.alignments.ConcatSortedAlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;

import java.io.IOException;

/**
 * A realignment processor that simply returns entries from the input reader.
 *
 * @author Fabien Campagne
 *         Date: May 12, 2011
 *         Time: 9:17:07 AM
 */
public class DummyProcessor implements AlignmentProcessorInterface {
    private ConcatSortedAlignmentReader reader;

    /**
     * Set the input reader.
     *
     * @param sortedReaders Input reader.
     */
    public DummyProcessor(ConcatSortedAlignmentReader sortedReaders) {
        this.reader = sortedReaders;
    }

    /**
     * Return the next available entry in the input reader, at targetIndex and position.
     *
     * @param targetIndex
     * @param position
     * @return
     * @throws IOException
     */
    public Alignments.AlignmentEntry nextRealignedEntry(int targetIndex, int position) throws IOException {
        return reader.skipTo(targetIndex, position);
    }

    public void setGenome(RandomAccessSequenceInterface genome) {

    }
}
