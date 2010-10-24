/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.util;

import it.unimi.dsi.fastutil.objects.ObjectHeapPriorityQueue;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import edu.cornell.med.icb.goby.methylation.MethylationSimilarityMatch;


/**
 * This class is an adaptation of MG4J's DocumentScoreBoundedSizeQueue.
 *
 * @author Fabien Campagne Date: Oct 21 2010
 */
public class HitBoundedPriorityQueue {

    /**
     * The underlying queue.
     */
    protected final ObjectHeapPriorityQueue<MethylationSimilarityMatch> queue;
    /**
     * The maximum number of documents to be ranked.
     */
    protected final int maxSize;


    /**
     * Creates a new empty bounded-size queue with a given capacity and natural
     * order as comparator. <p/> <P>Documents with equal scores will be compared
     * using their document index.
     *
     * @param capacity the initial capacity of this queue.
     */
    public HitBoundedPriorityQueue(final int capacity) {
        super();
        maxSize = capacity;
        queue = new ObjectHeapPriorityQueue<MethylationSimilarityMatch>(capacity, MethylationSimilarityMatch.INCREASING_SCORE_COMPARATOR);

    }


    /**
     * Checks whether a document with given score would be enqueued. <p/> <p>If
     * this methods returns true, an immediately following call to {@link
     * #enqueue(int,float)} with the same score would cause the document to be
     * enqueued.
     *
     * @param targetPosition position on the target
     * @param score          its score.
     * @return true if the match would have been enqueued by {@link
     *         #enqueue(int,float)}.
     */

    public boolean wouldEnqueue(final int targetPosition, final int chromosome, final float score) {
        if (queue.size() < maxSize && !targetPositions.contains(targetPosition)) {
            return true;
        }
        final MethylationSimilarityMatch dsi = queue.first();

        return score > dsi.score;
    }

    IntOpenHashSet targetPositions = new IntOpenHashSet();

    /**
     * Enqueues a transcript with given score and info.
     *
     * @param targetPosition position on the target
     * @param score          its score
     * @return true if the document has been actually enqueued.
     */

    public boolean enqueue(final int chromosome, final int targetPosition, final float score,
                           final int startForward,final int  endForward,final int  startReverse,
                           final int  endReverse,
                           final int windowLength,final float sumForwardStrand,
                           float sumReverseStrand) {

        if (maxSize == 0 || targetPositions.contains(targetPosition)) {
            return false;
        }
        if (queue.size() < maxSize) {
            if (targetPositions.contains(targetPosition)) {
                return false;
            }
            queue.enqueue(new MethylationSimilarityMatch(score, chromosome,  targetPosition));
            return true;
        } else {
            final MethylationSimilarityMatch dsi = queue.first();

            if (score > dsi.score) {
                targetPositions.remove(dsi.targetPosition);
                dsi.targetPosition = targetPosition;
                dsi.chromosome=chromosome;
                dsi.score = score;
                dsi.windowLength=windowLength;
                dsi.sumForwardStrand=sumForwardStrand;
                dsi.sumReverseStrand=sumReverseStrand;
                dsi.startForward = startForward;
                dsi.endForward = endForward;
                dsi.startReverse = startReverse;
                dsi.endReverse = endReverse;
                queue.changed();
                targetPositions.add(targetPosition);

                return true;
            }
            return false;
        }
    }
    public boolean isEmpty() {
        return queue.isEmpty();
    }

    public int size() {
        return queue.size();
    }

    /**
     * Dequeues a document from the queue, returning an instance of {@link
     * edu.cornell.med.icb.tissueinfo.similarity.TranscriptScore}.
     * Documents are dequeued in inverse score order.
     *
     * @return the next {@link edu.cornell.med.icb.tissueinfo.similarity.TranscriptScore}.
     */
    public final MethylationSimilarityMatch dequeue() {

        final MethylationSimilarityMatch transcriptScore = queue.dequeue();

        return transcriptScore;
    }


}
