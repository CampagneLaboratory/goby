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

import edu.cornell.med.icb.goby.alignments.ConcatSortedAlignmentReader;

/**
 * This factory returns the DummyProcessor.
 * @author Fabien Campagne
 *         Date: 5/31/11
 *         Time: 1:21 PM
 */
public class DefaultAlignmentProcessorFactory implements AlignmentProcessorFactory {
    /**
     * Always returns a DummyProcessor.
     * @param sortedReaders Input to the processor
     * @return a DummyProcessor that will not change the alignment.
     */
    public AlignmentProcessorInterface create(ConcatSortedAlignmentReader sortedReaders) {
        return new DummyProcessor(sortedReaders);
    }
}
