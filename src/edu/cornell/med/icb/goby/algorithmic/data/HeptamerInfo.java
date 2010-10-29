/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.algorithmic.data;

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.shorts.Short2FloatMap;
import it.unimi.dsi.fastutil.shorts.Short2FloatOpenHashMap;
import org.apache.commons.io.IOUtils;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * @author Fabien Campagne
 *         Date: May 17, 2010
 *         Time: 11:24:48 AM
 */
public class HeptamerInfo implements Serializable {
    private static final long serialVersionUID = -6209837092878262426L;

    public IndexedIdentifier heptamerToIndices = new IndexedIdentifier();
    public Short2FloatMap heptamerIndexToWeight = new Short2FloatOpenHashMap();
    public int heptamerLength = 7;
    public boolean colorSpace = false;
    /**
     * Load heptamer info from disk.
     *
     * @param filename
     * @return
     * @throws IOException
     * @throws ClassNotFoundException
     */
    public static HeptamerInfo load(final String filename) throws IOException, ClassNotFoundException {
        GZIPInputStream inputStream = null;
        try {
            inputStream = new GZIPInputStream(new FileInputStream(filename));
            return (HeptamerInfo) BinIO.loadObject(inputStream);
        } finally {
            if (inputStream != null) {
                IOUtils.closeQuietly(inputStream);

            }
        }
    }

    /**
     * Save heptamer info to disk.
     *
     * @param filename
     * @throws IOException
     * @throws ClassNotFoundException
     */

    public void save(final String filename) throws IOException {
        GZIPOutputStream gzipOutputStream = null;

        try {
            gzipOutputStream = new GZIPOutputStream(new FileOutputStream(filename));
            BinIO.storeObject(this, gzipOutputStream);
        } finally {
            IOUtils.closeQuietly(gzipOutputStream);
        }
    }
}
