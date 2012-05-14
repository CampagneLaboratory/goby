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
import edu.cornell.med.icb.goby.reads.ReadProtobuffCollectionHandler;
import edu.cornell.med.icb.goby.reads.Reads;
import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import static junit.framework.Assert.assertNotNull;

/**
 * @author Fabien Campagne
 *         Date: 5/12/12
 *         Time: 11:11 AM
 */
public class TestBZip2ChunkCodec {


    @Test
    public void roundTrip() throws IOException {
        BZip2ChunkCodec codec=new BZip2ChunkCodec();
        codec.setHandler(new ReadProtobuffCollectionHandler());
        Message collection= Reads.ReadCollection.newBuilder().build();
        ByteArrayOutputStream stream = codec.encode(collection);
        Message decoded = codec.decode(stream.toByteArray());
        assertNotNull(decoded);
    }
}
