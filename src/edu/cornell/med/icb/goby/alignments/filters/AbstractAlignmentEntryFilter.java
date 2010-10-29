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

package edu.cornell.med.icb.goby.alignments.filters;

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.IndexedIdentifier;

/**
 * Abstract class for AlignmentEntryFilter. These assist with merging alignment files.
 *
 * @author Fabien Campagne
 * @author Kevin Dorff
 */
public abstract class AbstractAlignmentEntryFilter {
    /**
     * Give the filter access to targets of the merged alignment.
     *
     * @param targets
     */
    public abstract void setTargetIdentifiers(final IndexedIdentifier targets);

    /**
     * Called during first pass processing of every entry.
     *
     * @param entry the entry to inspect for candidacy in the merge
     */
    public abstract void inspectEntry(final Alignments.AlignmentEntry entry);

    /**
     * Called after first pass (inspecEntry) processing is complete.
     */
    public abstract void postProcessing();

    /**
     * Returns true if the entry should be retained.
     * @param entry the entry to be inspected.
     * @return whether or not the entry should be retained
     */
    public abstract boolean shouldRetainEntry(final Alignments.AlignmentEntry entry);


    public void printStats() {

    }
}
