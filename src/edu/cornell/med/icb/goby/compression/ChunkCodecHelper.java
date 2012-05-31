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

import it.unimi.dsi.fastutil.bytes.ByteArraySet;
import it.unimi.dsi.fastutil.bytes.ByteSet;
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
                LOG.warn("Could not find codec " + codecName);
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

    /**
     * This method returns a codec associates with the registration code, or null when none is found. It does not raise an
     * error when the code does not match a known codec.
     *
     * @param registrationCode
     * @return
     */
    public static synchronized ChunkCodec withRegistrationCodeSilent(final byte registrationCode) {

        codecLoader.reload();
        for (final ChunkCodec chunkCodec : codecLoader) {
            if (chunkCodec.registrationCode() == registrationCode) {
                return chunkCodec;

            }
        }
        return null;

    }

    /**
     * Return the set of registration codes supported by this implementation.
     * @return the set of registration codes supported by this implementation.
     */
    public static synchronized ByteSet registrationCodes() {
        final ByteSet result=new ByteArraySet();
        codecLoader.reload();
        for (final ChunkCodec chunkCodec : codecLoader) {
            result.add(chunkCodec.registrationCode());
        }
        return result;
    }
}
