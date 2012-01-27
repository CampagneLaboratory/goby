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

import it.unimi.dsi.fastutil.ints.IntAVLTreeSet;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.algorithmic.data.UnboundedFifoPool;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;

import java.util.Random;

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
    public IntSet positionsWithSpanningIndel = new IntAVLTreeSet();
    UnboundedFifoPool<Alignments.AlignmentEntry> entriesInWindow = new UnboundedFifoPool<Alignments.AlignmentEntry>();
    public ObjectAVLTreeSet<ObservedIndel> potentialIndels = new ObjectAVLTreeSet<ObservedIndel>();
    /**
     * Position of the start of the realignment window for this target sequence. Alignment entries located
     * between windowStartPosition and windowStartPosition+windowLength are stored in entriesInWindow
     */
    public int windowStartPosition = Integer.MAX_VALUE;
    /**
     * The maximum entry position reached on the target sequence. If when maxEntryPosition >windowStartPosition+windowLength,
     * we do not need to keep adding to the pool, since there is enough in the pool to return another entry.
     */
    public int maxEntryPosition;
    /**
     * The number of entries encountered pas the maximum threshold.
     */
    protected int pastMaxCount;
    private static final int MAX_ENTRIES_IN_WINDOW = 500000;


    public void addIndel(int startPosition, int endPosition, String from, String to) {
        for (int p = startPosition; p < endPosition; p++) {
            positionsWithSpanningIndel.add(p);
        }
        final ObservedIndel candidate = new ObservedIndel(startPosition, endPosition, from, to);
        if (!potentialIndels.contains(candidate)) {
            potentialIndels.add(candidate);
            //    System.out.printf("Adding indel %s %n", candidate);
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

        for (final int pos : positionsWithSpanningIndel) {
            //  if (pos >= firstPosition && pos < lastPosition) {
            if (pos < lastPosition && pos < firstPosition) {

                positionsWithSpanningIndel.remove(pos);
            }
        }

        for (final ObservedIndel indel : potentialIndels) {

            if (indel.getEnd() <= lastPosition) {
                potentialIndels.remove(indel);

            }
        }
    }

    public void clear() {
        potentialIndels.clear();
        positionsWithSpanningIndel.clear();
        entriesInWindow.clear();
        pastMaxCount = 0;
    }

    public Alignments.AlignmentEntry remove() {
        final Alignments.AlignmentEntry tmp = entriesInWindow.remove();
        if (entriesInWindow.size() < MAX_ENTRIES_IN_WINDOW) {
            // when we have removed enough entries, make sure the next very high peak get the initial behavior:
            pastMaxCount = 0;
        }
        return tmp;
    }

    private static final long SPECIAL_SEED = 238927383682638267L;
    private Random random = new Random(SPECIAL_SEED);


    /**
     * Consider adding an entry. Entries are added unless we have reached the max count for the window.
     *
     * @param entry
     * @return True when the entry was added. False otherwise.
     */
    public boolean add(Alignments.AlignmentEntry entry) {
        boolean add = true;
        if (entriesInWindow.size() > MAX_ENTRIES_IN_WINDOW) {
            pastMaxCount++;
            final double randomChoice = random.nextDouble();
            final double threshold = 1.0 / pastMaxCount;
            if (randomChoice > threshold) {
                // we make it increasingly hard to add new entries past the max_entries threshold. This prevents from
                // running out of memory in the realignment step at positions that have clonal peaks (see this with RRBS).
                add = false;
            }
        }
        if (add) {
            entriesInWindow.add(entry);

        }
        return add;
    }
}
