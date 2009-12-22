/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.algorithmic.data;

import java.util.Comparator;

public class Read {
    public final int start;
    public final int end;

    public Read(final int start, final int end) {
        this.start = start;
        this.end = end;
    }

    public static final class ReadSortByStart implements Comparator<Read> {
        public int compare(final Read o1, final Read o2) {
            return o1.start - o2.start;
        }
    }

    public static final class ReadSortByEnd implements Comparator<Read> {
        public int compare(final Read o1, final Read o2) {
            return o1.end - o2.end;
        }
    }
}
