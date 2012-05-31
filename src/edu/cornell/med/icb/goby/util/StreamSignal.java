/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

import it.unimi.dsi.lang.MutableString;

import java.io.IOException;
import java.io.OutputStream;

/**
 *  A stream that checks if a token is written on a line and sets a flag when the event is observed.
 *  @author Fabien Campagne
 *         Date: 11/19/11
 *         Time: 6:27 PM
 */
public class StreamSignal extends OutputStream {
    MutableString scanning;
    private IsDone flag;
    private OutputStream stream;

    public StreamSignal(IsDone flag, String scanning, OutputStream stream) {
        this.scanning = new MutableString(scanning).compact();
        this.flag = flag;
        this.stream = stream;
    }

    /**
     * Store data until end of line.
     */

    MutableString string = new MutableString();

    @Override
    public void write(int b) throws IOException {
        if (b == '\n') {
            if (string.indexOf(scanning) >= 0) {
                flag.setDone(true);
            }
            string.setLength(0);

        }
        string.append((char) b);
        stream.write(b);

    }
}
