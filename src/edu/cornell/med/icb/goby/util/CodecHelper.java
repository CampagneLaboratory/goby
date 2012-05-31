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

package edu.cornell.med.icb.goby.util;

import com.google.protobuf.ByteString;
import edu.cornell.med.icb.goby.reads.ReadCodec;
import it.unimi.dsi.fastutil.bytes.Byte2ObjectMap;
import it.unimi.dsi.fastutil.bytes.Byte2ObjectOpenHashMap;

import java.util.ServiceLoader;

/**
 * @author Fabien Campagne
 *         Date: 1/11/12
 *         Time: 3:57 PM
 */
public class CodecHelper {
    private static final ServiceLoader<ReadCodec> readCodecLoader = ServiceLoader.load(ReadCodec.class);

    private static Byte2ObjectMap<ReadCodec> codeToReadCodec = new Byte2ObjectOpenHashMap(5);

    public static void reload() {
        readCodecLoader.reload();
        for (final ReadCodec next : readCodecLoader) {
            codeToReadCodec.put(next.registrationCode(), next);
        }
    }

    static {
        reload();
    }
    public static  ReadCodec locateReadCodec(final ByteString compressedData) {
        if (codeToReadCodec.size()==0) {
            reload();
        }
        final byte codecId = compressedData.byteAt(0);
        return codeToReadCodec.get(codecId);
    }
}
