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

package edu.cornell.med.icb.goby.reads;

import com.google.protobuf.GeneratedMessage;
import com.google.protobuf.Message;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * A mechanism to obtain a protobuff collection from a stream of data. Parsers must know the type of protobuff messages
 * contained in the collection.
 *
 * @author Fabien Campagne
 *         Date: 3/3/12
 *         Time: 11:37 AM
 */
public interface ProtobuffCollectionHandler {
    public final int TYPE_READS = 0;
    public final int TYPE_ALIGNMENTS = 1;

    /**
     * Returns the type of the collection elements.
     * @return One of the pre-defined types, TYPE_READS or TYPE_ALIGNMENTS.
     */
    public int getType();

    public GeneratedMessage parse(InputStream uncompressedStream) throws IOException;

    Message compressCollection(Message readCollection, ByteArrayOutputStream compressedBits) throws IOException;

    Message decompressCollection(Message reducedProtoBuff, byte[] compressedBytes);
}
