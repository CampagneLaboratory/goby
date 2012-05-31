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

import org.junit.Test;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;

import static org.junit.Assert.assertTrue;

/**
 * @author Fabien Campagne
 *         Date: 5/10/12
 *         Time: 1:47 PM
 */
public class TestChunkCodecs {


    @Test
    public void testValidateHybrid() throws IOException {

        DataInputStream dis = new DataInputStream(new FileInputStream("test-data/alignment-hybrid-codec/EJOYQAZ-small-hybrid.entries"));
        dis.skip(9);

        HybridChunkCodec1 hybridCodec = new HybridChunkCodec1();
        assertTrue(hybridCodec.validate((byte) 0, dis));

    }


    @Test
    public void testValidateGZip() throws IOException {

        DataInputStream dis = new DataInputStream(new FileInputStream("test-data/alignment-hybrid-codec/EJOYQAZ-small-gzip.entries"));
        dis.skip(9);

        GZipChunkCodec chunkCodec = new GZipChunkCodec();
        assertTrue(chunkCodec.validate((byte) 0, dis));

    }

    @Test
    public void testValidateBZip2() throws IOException {

        DataInputStream dis = new DataInputStream(new FileInputStream("test-data/alignment-hybrid-codec/EJOYQAZ-small-bzip2.entries"));
        dis.skip(9);

        BZip2ChunkCodec chunkCodec = new BZip2ChunkCodec();
        assertTrue(chunkCodec.validate((byte) 0, dis));

    }

}
