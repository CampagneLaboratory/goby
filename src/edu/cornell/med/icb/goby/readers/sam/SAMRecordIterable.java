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
 * If you use this in a foreach loop until the end of the loop (so that hasNext will be called until it
 * returns false) this will automatically close the iterator. If you need to exit early, you should call close().
 */
public class SAMRecordIterable implements CloseableIterator<SAMRecord>, Iterable<SAMRecord> {

    private final SAMRecordIterator base;

    private boolean closed;

    public SAMRecordIterable(final SAMRecordIterator base) {
        this.base = base;
        closed = false;
    }

    @Override
    public Iterator<SAMRecord> iterator() {
        return base;
    }

    /**
     * Manually (or automatically) close the base iterator.
     */
    @Override
    public void close() {
        if (!closed) {
            base.close();
        }
        closed = true;
    }

    /**
     * Determine if there is a next SAMRecord. If false, this will close the base iterator.
     * @return
     */
    @Override
    public boolean hasNext() {
        final boolean hasNext = base.hasNext();
        if (!hasNext) {
            close();
        }
        return hasNext;
    }

    /**
     * If you call next and this is already closed, it will return null.
     * @return the SAMRecord.
     */
    @Override
    public SAMRecord next() {
        if (closed) {
            return null;
        } else {
            return base.next();
        }
    }

    @Override
    public void remove() {
    }
}
