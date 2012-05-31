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

import it.unimi.dsi.fastutil.objects.ObjectListIterator;
import edu.cornell.med.icb.goby.alignments.Alignments;

/**
 * Provide a skipTo interface around an iterator over alignment entries. Useful for testing.
 *
 * @author Fabien Campagne
 *         Date: May 1, 2011
 *         Time: 12:13:16 PM
 */
public class SkipToListIterator extends SkipToIterator {
    private ObjectListIterator<Alignments.AlignmentEntry> iterator;
    /**
     * Initialize with the iterator.
     * @param iterator input iterator
     */
    public SkipToListIterator(ObjectListIterator<Alignments.AlignmentEntry> iterator) {
        this.iterator = iterator;
    }

    /**
     *
     * @param currentMinTargetIndex
     * @param position
     * @return
     */
    public Alignments.AlignmentEntry skipTo(int currentMinTargetIndex, int position) {
        while (iterator.hasNext()) {

            Alignments.AlignmentEntry element = iterator.next();
            if (element.getTargetIndex() >= currentMinTargetIndex && element.getPosition() >= position) return element;
        }
        return null;
    }
}
