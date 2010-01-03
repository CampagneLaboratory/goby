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

package edu.cornell.med.icb.nextgen.datamodel;

import java.io.Serializable;
import java.util.Date;
import java.util.Random;

/**
 * Uniquely tags a piece of information.
 *
 * @author Fabien Campagne
 *         Date: May 21, 2009
 *         Time: 11:49:18 AM
 */
public class Tagged implements Serializable {

    private static final long serialVersionUID = 154423600623605369L;

    public static int TAG_LENGTH_DEFAULT = 7;

    private String tag;

    private int intTag;

    private final int tagLength;

    public Tagged() {
        tagLength = TAG_LENGTH_DEFAULT;
    }

    public Tagged(final int tagLength) {
        assert tagLength >= 1;
        this.tagLength = tagLength;
    }

    /**
     * Get the tag. A unique identifier built by appending random uppercase letters up to a length.
     *
     * @return
     */
    public String getTag() {
        return tag;
    }

    /**
     * Get the tag. A unique integer value equivalent to the getTag() string representation. Equivalence means
     * that a.getTag().equals(b.getTag()) iff a.getIntTag() == b.getIntTag().
     *
     * @return
     */
    public int getIntTag() {
        return intTag;
    }

    static Random random = new Random(new Date().getTime());

    public void newTag() {
        final StringBuffer newTag = new StringBuffer();
        int newIntTag = 0;

        for (int i = 0; i < tagLength; i++) {
            //random character between A and Z:
            final int characterRange = (int) 'Z' - (int) 'A';
            final int randomValue = random.nextInt(characterRange) + 1;

            final char c = (char) ((int) 'A' + randomValue);
            newTag.append(c);
            final int shift=5;
            assert (Math.pow(2,shift))>characterRange :"character range is too large for bit shift value";
            newIntTag = (newIntTag << shift) | randomValue;
        }
        tag = newTag.toString();
        intTag = newIntTag;
    }

    /**
     * Useful when you only want a new String tag and don't need a full
     * Tagged object.
     * @return a new String tag
     */
    public static String createStringTag() {
        return createStringTag(TAG_LENGTH_DEFAULT);
    }

    public static String createStringTag(final int tagLength) {
        final Tagged tagged = new Tagged(tagLength);
        tagged.newTag();
        return tagged.getTag();
    }
}
