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

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;

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
    public static void main(String args[]) throws IOException {
        File directory = new File(".");

        String[] list = directory.list(new FilenameFilter() {
            public boolean accept(File directory, String filename) {
               
                final String extension = FilenameUtils.getExtension(filename);
                return (extension.equals("entries"));
            }
        });
        for (String filename : args) {
            String extension = FilenameUtils.getExtension(filename);
            String basename = FilenameUtils.removeExtension(filename);
            for (String alignFilename : list) {
                String alignBasename = FilenameUtils.removeExtension(alignFilename);
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
