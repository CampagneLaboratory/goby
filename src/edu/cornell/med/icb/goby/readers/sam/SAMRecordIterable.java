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

package edu.cornell.med.icb.goby.readers.sam;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloseableIterator;

import java.util.Iterator;

/**
 * Make a SAMRecordIterator iterable and thus easier to use.
 */
public class SAMRecordIterable implements CloseableIterator<SAMRecord>, Iterable<SAMRecord> {

    private final SAMRecordIterator base;

    public SAMRecordIterable(final SAMRecordIterator base) {
        this.base = base;
    }

    @Override
    public Iterator<SAMRecord> iterator() {
        return base;
    }

    @Override
    public void close() {
        base.close();
    }

    @Override
    public boolean hasNext() {
        return base.hasNext();
    }

    @Override
    public SAMRecord next() {
        return base.next();
    }

    @Override
    public void remove() {
    }
}
