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


import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.Writer;

/**
 * @author Fabien Campagne
 *         Date: 2/19/12
 *         Time: 4:38 PM
 */
public class OutputInfoFromWriter extends OutputInfo  {
    private Writer writer;
    public OutputInfoFromWriter(Writer writer) {
        this.writer=writer;

    }

    @Override
    public OutputStream getOutputStream() {
        throw new UnsupportedOperationException("Obtaining an output stream from this outputinfo implementation is not supported.");
    }

    @Override
    public PrintWriter getPrintWriter() {
        return new PrintWriter(writer);
    }

}
