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

package edu.cornell.med.icb.goby.alignments;

/**
 * Compression/Decompression (codec) to encode and decode compressed data in protocol buffer encoded alignments.
 *
 * @author Michael Meyer
 *         Date: 1/5/12
 *         Time: 2:19 PM
 */
public interface AlignmentCodec {
    /**
     * Encode an alignment into compressed form. The codec may compress only a subset of the fields of source. In this case,
     * the fields that cannot be compressed are left as protocol buffer fields and returned in the result alongside with
     * the compressed_data field.
     *
     * @param source The alignment to encode/compress.
     * @return The alignment with compressed data, if the codec could process source, null otherwise.
     */
    public Alignments.AlignmentEntry.Builder encode(Alignments.AlignmentEntry.Builder source);

    /**
     * Decode the compressed data of an alignment. Source is returned if the alignment does not contain any compressed data.
     *
     * @param source The alignment to decode (must contain the field compressed_data)
     * @return The decoded alignment, or null, if the codec could not handle the compressed data stored in source.
     */
    public Alignments.AlignmentEntry.Builder decode(Alignments.AlignmentEntry source);

    /**
     * Return the name of this codec.
     *
     * @return Return the name of this codec.
     */
    public String name();

    /**
     * This method is called when a new chunk of data is started. The event can be used by some codecs to reset
     * their probabilistic models of symbols.
     */
    public void newChunk();

    /**
     * Return the registration code of this codec, a byte that uniquely identifies this codec.
     *
     * @return Return the registration code of this codec.
     */
    byte registrationCode();
}
