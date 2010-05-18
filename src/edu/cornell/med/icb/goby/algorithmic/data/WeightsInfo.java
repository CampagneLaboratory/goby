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

import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.floats.FloatArrayList;

import java.io.IOException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.Serializable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.reads.ReadsReader;

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

    FloatArrayList weights = new FloatArrayList();

    public void setWeight(int readIndex, float weight) {
        if (weights.size() - 1 < readIndex) {
            weights.size((readIndex + 10) * 2);
        }

        weights.set(readIndex, weight);
    }

    public float getWeight(int readIndex) {
        return weights.getFloat(readIndex);
    }

    /**
     * Load weights info from disk.
     *
     * @param filename
     * @return
     * @throws java.io.IOException
     * @throws ClassNotFoundException
     */
    public static WeightsInfo load(String filename) throws IOException, ClassNotFoundException {
        GZIPInputStream inputStream = null;
        try {
            inputStream = new GZIPInputStream(new FileInputStream(filename));
            return (WeightsInfo) BinIO.loadObject(inputStream);
        } finally {
            if (inputStream != null) {
                IOUtils.closeQuietly(inputStream);

            }
        }
    }

    /**
     * Tries to load a weight information file matching basename. Simply takes out the file extension and
     * appends .hepmap to try and locate the weights corresponding to the reads or alignment file.
     *
     * @param filename name of a compact reads or alignment
     * @return
     * @throws IOException
     * @throws ClassNotFoundException
     */
    public static WeightsInfo loadForBasename(String filename) throws IOException, ClassNotFoundException {
        String basename = AlignmentReader.getBasename(filename);
        basename = ReadsReader.getBasename(basename);
        return load(basename + ".weights");
    }

    /**
     * Save weights info to disk.
     *
     * @param filename
     * @throws IOException
     * @throws ClassNotFoundException
     */

    public void save(String filename) throws IOException {

        GZIPOutputStream gzipOutputStream = null;

        try {
            gzipOutputStream = new GZIPOutputStream(new FileOutputStream(filename));
            BinIO.storeObject(this, gzipOutputStream);
            gzipOutputStream.flush();
            LOG.info("Saved " + filename);
        }
        finally {
            IOUtils.closeQuietly(gzipOutputStream);
        }

    }

    public void size(int numberOfReads) {
        weights.size(numberOfReads);
    }
}
