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

package edu.cornell.med.icb.goby.methylation;

import java.util.Comparator;

/**
 * @author Fabien Campagne
*         Date: Oct 22, 2010
*         Time: 12:13:47 PM
*/
public class MethylationSimilarityMatch {
    public float score;
    public int targetPosition;
    public int chromosome;
    public int windowLength;
    public float sumReverseStrand;
    public float sumForwardStrand;
    public int startForward;
    public int endForward;
    public int startReverse;
    public int endReverse;

    public MethylationSimilarityMatch(float score, int chromosome, int targetPosition) {
        this.score = score;
        this.chromosome=chromosome;
        this.targetPosition = targetPosition;
    }

    public static Comparator<MethylationSimilarityMatch> INCREASING_SCORE_COMPARATOR = new Comparator<MethylationSimilarityMatch>() {
        public int compare(MethylationSimilarityMatch a, MethylationSimilarityMatch b) {
            return Float.compare(a.score, b.score);
        }
    };
}
