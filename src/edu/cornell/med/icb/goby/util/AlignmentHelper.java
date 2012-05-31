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

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;

import java.io.File;

/**
 * @author Fabien Campagne
 *         Date: 4/17/12
 *         Time: 6:45 PM
 */
public class AlignmentHelper {

    /**
     * Collect the set of files that make up the alignment output of the aligner.
     *
     * @param basename Basename where the alignment was written.
     * @return an array of files that will contain the output from the aligner
     */
    public static File[] buildResults(final String basename) {
        final ObjectArrayList<File> result = new ObjectArrayList<File>();
        final String[] extensions = {".stats", ".entries", ".header", ".tmh"};

        for (final String extension : extensions) {
            final File file = new File(basename + extension);
            if (file.exists()) {
                result.add(file);
            }
        }

        return result.toArray(new File[result.size()]);
    }


}
