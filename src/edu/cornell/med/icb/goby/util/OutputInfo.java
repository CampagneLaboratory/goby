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

import java.io.*;

/**
 * @author Fabien Campagne
 *         Date: 2/19/12
 *         Time: 4:38 PM
 */
public class OutputInfo {
    private PrintWriter outWriter;
    private OutputStream outStream;
    private String outputFilename;


    public OutputInfo(String outputFilename) throws FileNotFoundException {
        this.outputFilename = outputFilename;
    }

    public OutputInfo() {
    }

    public PrintWriter getPrintWriter() {
        try {
            outWriter = isToConsole(outputFilename) ? new PrintWriter(System.out) : new PrintWriter(outputFilename);
        } catch (FileNotFoundException e) {
            System.err.println("Cannot open output file for writing: " + outputFilename);
        }
        return outWriter;
    }

    public Writer getWriter() {
        return null;
    }

    public boolean isToConsole(String outputFilename) {
        return "-".equals(outputFilename);
    }

    public OutputStream getOutputStream() {
        try {
            outStream = new FileOutputStream(outputFilename);
        } catch (FileNotFoundException e) {
            System.err.println("Cannot open output file for writing: " + outputFilename);
        }
        return outStream;
    }

    public String getFilename() {
        return outputFilename;
    }
}
