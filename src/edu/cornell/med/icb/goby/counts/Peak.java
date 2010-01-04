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

package edu.cornell.med.icb.goby.counts;

/**
 * Information about a peak of something (e.g., read counts) along a reference.
 *
 * @author Fabien Campagne
 *         Date: May 27, 2009
 *         Time: 6:44:06 PM
 */
public class Peak {
    public int start;
    public int count;
    public int length;

    /**
     * Create a copy of this peak information.
     *
     * @return a new instance of Peak with the same information.
     */
    public Peak copy() {
        final Peak result = new Peak();
        result.start = this.start;
        result.count = this.count;
        result.length = this.length;
        return result;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append(" peak :");
        sb.append(" start :");
        sb.append(start);
        sb.append(" count :");
        sb.append(count);
        sb.append(" length :");
        sb.append(length);
        return sb.toString();
    }
}
