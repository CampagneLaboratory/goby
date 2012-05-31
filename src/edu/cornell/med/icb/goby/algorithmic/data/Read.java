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

package edu.cornell.med.icb.goby.algorithmic.data;

import java.io.Serializable;
import java.util.Comparator;

public class Read {
    public final int start;
    public final int end;



    /**
     *
     * @param start Start position of the alignment entry.
     * @param end  End position of the alignment entry on the reference sequence.

     */
    public Read(final int start, final int end) {
        this.start = start;
        this.end = end;

    }

    public static final class ReadSortByStart implements Comparator<Read>, Serializable {
        private static final long serialVersionUID = 7938501857256287597L;

        public int compare(final Read o1, final Read o2) {
            return o1.start - o2.start;
        }
    }

    public static final class ReadSortByEnd implements Comparator<Read>, Serializable {
        private static final long serialVersionUID = 6267512969132795775L;

        public int compare(final Read o1, final Read o2) {
            return o1.end - o2.end;
        }
    }
}
