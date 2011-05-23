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

package edu.cornell.med.icb.goby.counts;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: May 23, 2011
 *         Time: 4:09:04 PM
 */
public class PositionAdapter implements CountsReaderI {
    /**
     * Return the position after the transition occured.
     * @return
     */
    public int getPosition() {
        return reader.getPosition()+reader.getLength();
    }

    public boolean hasNextTransition() throws IOException {
        return reader.hasNextTransition();
    }

    public void nextTransition() throws IOException {
        reader.nextTransition();
    }

    public int getCount() {
        return reader.getCount();
    }

    public void skipTo(int position) throws IOException {
        reader.skipTo(position);
    }

    public int getLength() {
        return reader.getLength();
    }

    public void close() throws IOException {
        reader.close();
    }

    public PositionAdapter(CountsReaderI reader) {
        this.reader = reader;
    }

    CountsReaderI reader;

}
