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


/**
 * @author Fabien Campagne
 *         Date: 6/15/11
 *         Time: 8:09 AM
 */
public abstract class AbstractConcatAlignmentReader extends AbstractAlignmentReader {
    /**
     * Whether this concat reader should record or modify the sampleIndex field of each entry according to the reader of
     * origin.
     */
    protected boolean adjustSampleIndices;

    /**
     * This constructor will upgrade the alignment to the latest version of the Goby data structure.
     *
     * @param upgrade  indicate whether the upgrade should be considered.
     * @param basename Alignment basename.
     */
    protected AbstractConcatAlignmentReader(final boolean upgrade, final String basename) {
        super(upgrade, basename);
    }



    /**
     * Instruct this concat reader to set sample indices on entries to reflect the reader of origin before concatenation.
     *
     * @param adjust When true, the concat will set sample indices, when false, sample indices are not modified or created.
     */
    public void setAdjustSampleIndices(final boolean adjust) {
        adjustSampleIndices = adjust;
    }

}
