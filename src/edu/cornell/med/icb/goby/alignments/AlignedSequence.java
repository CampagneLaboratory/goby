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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.lang.MutableString;

/**
 * Encodes information about one aligned sequence (support for MAF alignment format).
 *
 * @author Fabien Campagne
 *         Date: May 4, 2009
 *         Time: 3:54:17 PM
 */
public class AlignedSequence {
    public AlignedSequence() {
        sequenceIdentifier = new MutableString();
        alignment = new MutableString();
    }

    public MutableString sequenceIdentifier;
    public int alignedStart;
    public int alignedLength;
    public char strand;
    public int sequenceLength;
    public MutableString alignment;

}
