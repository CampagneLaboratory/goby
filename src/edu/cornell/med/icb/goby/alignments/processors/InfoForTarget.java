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

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.algorithmic.data.UnboundedFifoPool;

/**
 * @author Fabien Campagne
 *         Date: May 14, 2011
 *         Time: 5:59:16 PM
 */
public class InfoForTarget {
    int targetIndex;
    /**
     * positionsWithSpanningIndel contains zero-based positions.
     */
    public IntArraySet positionsWithSpanningIndel = new IntArraySet();
    UnboundedFifoPool<Alignments.AlignmentEntry> entriesInWindow = new UnboundedFifoPool<Alignments.AlignmentEntry>();
    public ObjectArrayList<ObservedIndel> potentialIndels = new ObjectArrayList<ObservedIndel>();
    /**
     * Position of the start of the realignment window for this target sequence. Alignment entries located
     * between windowStartPosition and windowStartPosition+windowLength are stored in entriesInWindow
     */
    public int windowStartPosition=Integer.MAX_VALUE;

    public void addIndel(int startPosition, int endPosition, String from, String to) {
        for (int p = startPosition; p < endPosition; p++) {
            positionsWithSpanningIndel.add(p);
        }
        final ObservedIndel candidate = new ObservedIndel(startPosition, endPosition, from, to);
        if (!potentialIndels.contains(candidate)) {
            potentialIndels.add(candidate);
         System.out.printf("Adding indel %s %n", candidate);
        }

    }

    public InfoForTarget(int targetIndex) {
        this.targetIndex = targetIndex;
    }


    /**
     * Remove indels that span firstPosition-lastPosition or whose end occur before endPosition (inclusive).
     *
     * @param firstPosition
     * @param lastPosition
     */
    public void removeIndels(int firstPosition, int lastPosition) {

        for (int pos : positionsWithSpanningIndel) {
            if (pos >= firstPosition && pos < lastPosition) {

               positionsWithSpanningIndel.remove(pos);
            }
        }


        // TODO is this removing enough?
        final ObjectArrayList<ObservedIndel> toRemove = new ObjectArrayList<ObservedIndel>();
        for (ObservedIndel indel : potentialIndels) {
            if (indel.getEnd() <= lastPosition) {

               toRemove.add(indel);
                System.out.printf("lastPosition=%d removing indel %s %n",lastPosition, indel);
            }
        }
        potentialIndels.removeAll(toRemove);

    }

    public void clear() {
        potentialIndels.clear();
        positionsWithSpanningIndel.clear();
        entriesInWindow.clear();
    }
}
