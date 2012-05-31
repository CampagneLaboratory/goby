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

package edu.cornell.med.icb.goby.reads;

import it.unimi.dsi.lang.MutableString;

import java.io.IOException;

/**
 * Provides access to either a Goby random access genome, or a Picard indexed fasta file.
 *
 * @author Fabien Campagne
 *         Date: 4/1/12
 *         Time: 3:00 PM
 */
public class DualRandomAccessSequenceCache implements RandomAccessSequenceInterface {
    /**
     * Load either a Picard fasta index genome, or goby random access genome. To load a picard
     * genome, you must provide a filename ending in .fa or .fasta. The fasta file must be indexed
     * with samtools faidx and the fasta index file must be named filename.fasta.idx or filename.fa.idx
     * To load a Goby random access cache, you need to create the cache with the goby build-sequence-cache
     * mode and provide a basename as filename.
     *
     * @param filename
     * @throws IOException
     */
    public void load(final String filename) throws IOException, ClassNotFoundException {
        if (filename.endsWith(".fa") || filename.endsWith(".fasta")) {

            delegate = new PicardFastaIndexedSequence(filename);
        } else {

            final RandomAccessSequenceCache gobyCache = new RandomAccessSequenceCache();
            gobyCache.load(filename);
            delegate = gobyCache;

        }

    }

    @Override
    public char get(int referenceIndex, int position) {
        return delegate.get(referenceIndex, position);
    }

    @Override
    public int getLength(int targetIndex) {
        return delegate.getLength(targetIndex);
    }

    @Override
    public void getRange(int referenceIndex, int position, int length, MutableString bases) {
        delegate.getRange(referenceIndex, position, length, bases);
    }

    @Override
    public int getReferenceIndex(String referenceId) {
        return delegate.getReferenceIndex(referenceId);
    }

    @Override
    public String getReferenceName(int index) {
        return delegate.getReferenceName(index);
    }

    @Override
    public int size() {
        return delegate.size();
    }

    private RandomAccessSequenceInterface delegate;


}
