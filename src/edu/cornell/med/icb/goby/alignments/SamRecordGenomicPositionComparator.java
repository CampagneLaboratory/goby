/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

import net.sf.samtools.SAMRecord;

import java.io.Serializable;
import java.util.Comparator;

/**
 * @author Fabien Campagne
 *         Date: 5/26/12
 *         Time: 4:53 PM
 */
public class SamRecordGenomicPositionComparator implements Comparator< SAMRecord>,Serializable {

    private static final long serialVersionUID = 3386079334589738180L;

    @Override
    public int compare(SAMRecord a,  SAMRecord b) {
       if (a.getReferenceIndex()==b.getReferenceIndex()) {
           return a.getAlignmentStart()-b.getAlignmentStart();
       }  else {
           return a.getReferenceIndex()-b.getReferenceIndex();
       }
    }
}
