/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.io.BinIO;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * @author Fabien Campagne
 *         Date: May 17, 2010
 *         Time: 1:10:19 PM
 */
public class WeightsInfo implements Serializable {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(WeightsInfo.class);
    private static final long serialVersionUID = 3965165699293207843l;

    protected final FloatArrayList weights = new FloatArrayList();

    public void setWeight(final int readIndex, final float weight) {
        if (weights.size() - 1 < readIndex) {
            weights.size((readIndex + 10) * 2);
        }

        weights.set(readIndex, weight);
    }

    public float getWeight(final int readIndex) {
        return weights.getFloat(readIndex);
    }

    /**
     * Load weights info from disk.
     *
     * @param filename The name of the file to load
     * @return The populated WeightsInfo object read from the file
     * @throws IOException If the file cannot be read
     * @throws ClassNotFoundException if the file contains a class that cannot be found.
     */
    public static WeightsInfo load(final String filename) throws IOException, ClassNotFoundException {
        GZIPInputStream inputStream = null;
        try {
            inputStream = new GZIPInputStream(new FileInputStream(filename));
            return (WeightsInfo) BinIO.loadObject(inputStream);
        } finally {
            IOUtils.closeQuietly(inputStream);
        }
    }

    /**
     * Tries to load a weight information file matching basename. The file is the basename with
     * the extension in the form of the id plus the string "-weights".
     *
     * @param filename name of a compact reads or alignment
     * @param id The id for the weights file to load
     * @return The populated WeightsInfo object read from the file
     * @throws IOException If the file cannot be read
     * @throws ClassNotFoundException if the file contains a class that cannot be found.
     */
    public static WeightsInfo loadForBasename(final String filename, final String id)
            throws IOException, ClassNotFoundException {
        // strip any compact alignment extensions
        String basename = AlignmentReader.getBasename(filename);

        // strip any compact reads extensions
        basename = ReadsReader.getBasename(basename);
        return load(basename + "." + id + "-weights");
    }

    /**
     * Save weights info to disk.
     *
     * @param filename The name of the file to save to
     * @throws IOException If the file cannot be written
     */
    public void save(final String filename) throws IOException {
        GZIPOutputStream gzipOutputStream = null;
        try {
            gzipOutputStream = new GZIPOutputStream(new FileOutputStream(filename));
            BinIO.storeObject(this, gzipOutputStream);
            gzipOutputStream.flush();
            LOG.info("Saved " + filename);
        } finally {
            IOUtils.closeQuietly(gzipOutputStream);
        }
    }

    public void size(final int numberOfReads) {
        weights.size(numberOfReads);
    }

    /**
     * Returns the number of weights stored.
     * @return maximum read index stored.
     */
    public int size() {
        return weights.size();
    }
}
