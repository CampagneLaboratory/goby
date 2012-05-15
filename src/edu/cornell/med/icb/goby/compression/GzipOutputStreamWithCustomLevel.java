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

import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.GZIPOutputStream;

/**
 * A sub-class of GZipOutputStream that makes it possible to configure the compression level.
 *
 * @author Fabien Campagne
 *         Date: 5/15/12
 *         Time: 2:45 PM
 */
public final class GzipOutputStreamWithCustomLevel extends GZIPOutputStream {

    public GzipOutputStreamWithCustomLevel(int level, final OutputStream out) throws IOException {
        super(out);
        setLevel(level);
    }

    /**
     * Set the level of GZip compression required.
     *
     * @param level A number between 0 and 9.
     */
    public void setLevel(int level) {
        def.setLevel(level);
    }
}