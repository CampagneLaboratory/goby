/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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
    /**
     * The name of one of the source sequences for the alignment.
     */
    public final MutableString sequenceIdentifier;

    /**
     * The start of the aligning region in the source sequence. This is a zero-based number.
     * If the @link #strand} is '-' then this is the start relative to the reverse-complemented
     * source sequence.
     */
    public int alignedStart;

    /**
     * The start of the aligning region in the source sequence. This is a zero-based number.
     * If the {@link #strand} is '-' then this is the start relative to the
     * reverse-complemented source sequence.
     */
    public int alignedLength;

    /**
     * Either '+' or '-'. If '-', then the alignment is to the reverse-complemented source.
     */
    public char strand;

    /**
     * The size of the entire source sequence, not just the parts involved in the alignment.
     */
    public int sequenceLength;

    /**
     * The nucleotides (or amino acids) in the alignment and any insertions (dashes) as well.
     */
    public final MutableString alignment;

    /**
     * Create a new aligned sequence object.
     * Typically contents of this object maps to Lines starting with 's'
     * (a sequence within an alignment block) in MAF alignment file as in the following example:
     *
     * <pre>
     * s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
     * s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
     * s baboon         249182 13 +   4622798 gcagctgaaaaca
     * s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
     * </pre>
     */
    public AlignedSequence() {
        super();
        sequenceIdentifier = new MutableString();
        alignment = new MutableString();
    }
}
