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

import com.google.protobuf.CodedInputStream;
import com.google.protobuf.GeneratedMessage;
import com.google.protobuf.Message;
import edu.cornell.med.icb.goby.compression.ProtobuffCollectionHandler;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * A handler for collections that contain read entries.
 *
 * @author Fabien Campagne
 *         Date: 3/3/12
 *         Time: 11:48 AM
 */
public class ReadProtobuffCollectionHandler implements ProtobuffCollectionHandler {
    @Override
    public int getType() {
        return TYPE_READS;
    }

    @Override
    public GeneratedMessage parse(final InputStream compressedBytes) throws IOException {
        final CodedInputStream codedInput = CodedInputStream.newInstance(compressedBytes);
        codedInput.setSizeLimit(Integer.MAX_VALUE);

        return Reads.ReadCollection.parseFrom(compressedBytes);
    }

    @Override
    public Message compressCollection(Message readCollection, ByteArrayOutputStream compressedBits) {
        return readCollection;
    }

    @Override
    public Message decompressCollection(Message reducedProtoBuff, byte[] compressedBytes) {
        return reducedProtoBuff;
    }

    @Override
    public void setUseTemplateCompression(boolean useTemplateCompression) {

    }
}
