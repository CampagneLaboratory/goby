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

package edu.cornell.med.icb.goby.maq;

/**
 * MAQ constants.
 *
 * @author Kevin Dorff
 */
public final class MaqConstants {

    /** Private constructor for static, constants class. */
    private MaqConstants() {
    }

    /** Default qual score for importing Eland data. */
    public static final short DEFAULT_QUAL = 25;

    /** Max name length. */
    public static final int MAX_NAMELEN = 36;

    // This is specified in the reader / writer classes.
    // public final static int MAX_READLEN = 64; // or 128

    /** MAQ MAP old format version number. */
    public static final int MAQMAP_FORMAT_OLD = 0;

    /** MAQ MAP new format version number. */
    public static final int MAQMAP_FORMAT_NEW = -1;

    /** Pair flag FF. */
    public static final int PAIRFLAG_FF = 0x01;

    /** Pair flag FR. */
    public static final int PAIRFLAG_FR = 0x02;

    /** Pair flag RF. */
    public static final int PAIRFLAG_RF = 0x04;

    /** Pair flag RR. */
    public static final int PAIRFLAG_RR = 0x08;

    /** Pair flag paired. */
    public static final int PAIRFLAG_PAIRED = 0x10;

    /** Pair flag diffchr. */
    public static final int PAIRFLAG_DIFFCHR = 0x20;

    /** Pair flag nowatch. */
    public static final int PAIRFLAG_NOMATCH = 0x40;

    /** Pair flag SW. */
    public static final int PAIRFLAG_SW = 0x80;

    /**
     * NST_NT4_TABLE for calculating Eland MAP information.
     */
    public static final short[] NST_NT4_TABLE = new short[] {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 };

    /**
     * Precalculate LOG_N values for speed.
     */
    public static final int[] LOG_N = new int[256];
    static {
        LOG_N[0] = -1;
        for (int i = 1; i != 256; ++i) {
            LOG_N[i] = (int) (3.434 * Math.log((float) i) + 0.5);
        }
    }
}
