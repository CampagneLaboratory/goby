/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

import it.unimi.dsi.fastutil.booleans.BooleanArrayList;
import it.unimi.dsi.fastutil.booleans.BooleanList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.logging.ProgressLogger;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Jun 20, 2009
 *         Time: 10:26:58 AM
 */
public class ReadsLoader {
    private final ReadSet readIndexFilter;
    private SequenceDigests[] digests;
    private final ObjectArrayList<byte[]> compressedReads;
    private int numReads;
    private final File readsFile;
    private static final int MAX_PROCESS_READS = Integer.MAX_VALUE;

    private int readLength;
    private byte[] byteBuffer;
    private final int numberOfMismaches = 0;

    private boolean colorSpace;

    public int getReadLength() {
        return readLength;
    }

    public ReadsLoader(final ReadSet readIndexFilter, final File readsFile) {
        super();
        this.readIndexFilter = readIndexFilter;
        compressedReads = new ObjectArrayList<byte[]>(1000000);
        this.readsFile = readsFile;
    }

    public void setColorSpace(final boolean colorSpace) {
        this.colorSpace = colorSpace;
    }

    public int read() throws IOException {
        final ProgressLogger progress = new ProgressLogger();
        progress.displayFreeMemory = true;
        progress.start("parsing reads");
        final ReadsReader readsReader = new ReadsReader(new FileInputStream(readsFile));

        int maxReadLength = 0;
        int minReadLength = Integer.MAX_VALUE;
        for (final Reads.ReadEntry readEntry : readsReader) {

            if (readIndexFilter == null || readIndexFilter.contains(readEntry.getReadIndex())) {
                maxReadLength = Math.max(maxReadLength, readEntry.getReadLength());
                minReadLength = Math.min(minReadLength, readEntry.getReadLength());
                if (minReadLength != maxReadLength) {
                    System.err.println("The read length must be fixed.");
                    System.exit(1);
                }
                if (digests == null) {
                    if (numberOfMismaches == 0) {
                        digests = new SequenceDigests[2];
                        digests[0] = new SequenceDigests(minReadLength, true);     // forward strand
                        digests[1] = new SequenceDigests(minReadLength, false);    // reverse strand

                    } else if (numberOfMismaches == 1) {
                        digests = new SequenceDigests[4];
                        final BooleanList mask1 = new BooleanArrayList();
                        final BooleanList mask2 = new BooleanArrayList();
                        for (int i = 0; i < readLength; i++) {

                            mask1.set(i, i % 2 == 1);
                            mask2.set(i, !mask1.get(i));
                        }
                        digests[0] = new SequenceDigests(mask1, minReadLength, true);
                        digests[1] = new SequenceDigests(mask1, minReadLength, false);
                        digests[2] = new SequenceDigests(mask2, minReadLength, true);
                        digests[3] = new SequenceDigests(mask2, minReadLength, false);
                    }
                }

                byteBuffer = new byte[maxReadLength]; //2 bits per base require four times less space.

                final int readIndex = readEntry.getReadIndex();
                final byte[] bytes = readEntry.getSequence().toByteArray();
                for (final SequenceDigests digest : digests) {

                    digest.digestAndStore(bytes, 0, readIndex);
                }
                System.arraycopy(bytes, 0, byteBuffer, 0, maxReadLength);

                if (readIndex > compressedReads.size() - 1) {
                    compressedReads.size((readIndex + 1) * 3 / 2);
                }
                compressedReads.set(readIndex, byteBuffer);

                progress.lightUpdate();

                if (numReads > MAX_PROCESS_READS) {
                    break;
                }
            }
            numReads++;
        }

        readLength = minReadLength;
        compressedReads.size(numReads + 1);
        readsReader.close();
        progress.stop("Finished parsing reads.");
        return numReads;
    }

    public ObjectList<byte[]> getCompressedReads() {
        return compressedReads;
    }

    public byte[] getByteBuffer() {
        return byteBuffer;
    }

    public SequenceDigests[] getDigests() {
        return digests;
    }
}
