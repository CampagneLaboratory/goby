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

package edu.cornell.med.icb.goby.util;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;

/**
 * Utility to rename weight files to matching alignment basenames.
 * Usage: java -cp goby.jar edu.cornell.med.icb.goby.util.RenameWeights *XXX.gc-weights
 * The previous line will copy each XXX/match/.gc-weight to the alignment basename
 * that matches XXX and ends with .entries.  For instance, the file XXX.gc-weights
 * will be renamed 1212-XXX.gc-weights if an alignment file named 1212-XXX.entries exists.
 *
 * @author Fabien Campagne
 *         Date: Jun 2, 2010
 *         Time: 11:34:32 AM
 */
public class RenameWeights {
    private RenameWeights() {
    }

    public static void main(final String[] args) throws IOException {
        final File directory = new File(".");

        final String[] list = directory.list(new FilenameFilter() {
            public boolean accept(final File directory, final String filename) {

                final String extension = FilenameUtils.getExtension(filename);
                return (extension.equals("entries"));
            }
        });
        for (final String filename : args) {
            final String extension = FilenameUtils.getExtension(filename);
            final String basename = FilenameUtils.removeExtension(filename);
            for (final String alignFilename : list) {
                final String alignBasename = FilenameUtils.removeExtension(alignFilename);
                if (alignBasename.endsWith(basename)) {
                    System.out.println("move " + filename + " to " + alignBasename + "." + extension);

                    final File destination = new File(alignBasename + "." + extension);
                    FileUtils.deleteQuietly(destination);
                    FileUtils.moveFile(new File(filename), destination);
                }
            }


        }
    }
}
