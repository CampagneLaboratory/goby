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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.util.ServiceLoader;

/**
 * A helper to load chunk codecs by name or registration code. The helper locates the codec as a
 * service. Services must be registered in the manifest of the jar file to be located.
 *
 * @author Fabien Campagne
 *         Date: 3/8/12
 *         Time: 8:35 AM
 */
public class ChunkCodecHelper {
    /**
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(ChunkCodecHelper.class);


    private static final ServiceLoader<ChunkCodec> codecLoader = ServiceLoader.load(ChunkCodec.class);

    private ChunkCodecHelper() {
    }

    public static synchronized ChunkCodec load(String codecName) {
        ChunkCodec codec = null;
        if (codecName != null) {
            codecLoader.reload();
            for (final ChunkCodec c : codecLoader) {

                if (c.name().equals(codecName)) {
                    LOG.debug("Will use chunk codec " + c.name());
                    codec = c;
                    break;
                }
            }
            if (codec == null) {
                LOG.info("Could not find codec " + codecName);
            }
        }
        return codec;
    }

    public static synchronized ChunkCodec withRegistrationCode(final byte registrationCode) {

        codecLoader.reload();
        for (final ChunkCodec chunkCodec : codecLoader) {
            if (chunkCodec.registrationCode() == registrationCode) {
                return chunkCodec;

            }
        }
        throw new InternalError("Codec registration code not recognized: " + registrationCode);

    }
}
