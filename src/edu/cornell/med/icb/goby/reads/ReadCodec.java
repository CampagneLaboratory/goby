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

package edu.cornell.med.icb.goby.reads;

/**
 * Compression/Decompression (codec) to encode and decode compressed data in protocol buffer encoded reads.
 *
 * @author Fabien Campagne
 *         Date: 8/19/11
 *         Time: 1:48 PM
 */
public interface ReadCodec {
    /**
     * Encode a read into compressed form. The codec may compress only a subset of the fields of source. In this case,
     * the fields that cannot be compressed are left as protocol buffer fields and returned in the result alongside with
     * the compressed_data field.
     *
     * @param source The read to encode/compress.
     * @return The read with compressed data, if the codec could process source, null otherwise.
     */
    public Reads.ReadEntry.Builder encode(Reads.ReadEntry.Builder source);

    /**
     * Decode the compressed data of a read. Source is returned if the read does not contain any compressed data.
     *
     * @param source The read to decode (must contain the field compressed_data)
     * @return The decoded read, or null, if the codec could not handle the compressed data stored in source.
     */
    public Reads.ReadEntry.Builder decode(Reads.ReadEntry source);

    /**
     * Return the name of this codec.
     * @return  Return the name of this codec.
     */
    public String name();

    /**
     * This method is called when a new chunk of data is started. The event can be used by some codecs to reset
     * their probabilistic models of symbols.
     */
    public void newChunk();

}
