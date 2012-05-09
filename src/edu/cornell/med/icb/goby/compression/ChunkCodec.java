/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.compression;

import com.google.protobuf.Message;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

/**
 * Encode/decode whole chunks of protocol buffer data to a stream of bytes.
 *
 * @author Fabien Campagne
 *         Date: 3/3/12
 *         Time: 10:27 AM
 */
public interface ChunkCodec {
    /**
     * Return the name of this codec.
     *
     * @return Return the name of this codec.
     */
    public String name();


    /**
     * Return the registration code of this codec, a byte that uniquely identifies this codec.
     *
     * @return Return the registration code of this codec.
     */
    byte registrationCode();

    /**
     * Encode the protobuff collection to a byte of stream.
     *
     * @param readCollection
     */
    ByteArrayOutputStream encode(Message readCollection) throws IOException;

    Message decode(byte[] bytes) throws IOException;

    /**
     * Return a suggestion for chunk size, to be used as default.
     * @return
     */
    public int getSuggestedChunkSize();

    public void setHandler(ProtobuffCollectionHandler parser) ;

    /**
     * Validate that the input stream is at the start of a chunk encoded with this codec.
     *
     * @param input input stream positioned at the start of the chunk.
     * @return True if the chunk is valid, false otherwise.
     */
    boolean validate(FastBufferedInputStream input);
}
