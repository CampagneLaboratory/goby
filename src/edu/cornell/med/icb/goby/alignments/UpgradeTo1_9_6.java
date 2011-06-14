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

package edu.cornell.med.icb.goby.alignments;

import com.google.protobuf.CodedInputStream;
import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FileUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * @author Fabien Campagne
 *         Date: May 6, 2011
 *         Time: 6:53:21 PM
 */
public class UpgradeTo1_9_6 {
    private boolean verbose;

    public void upgrade(String basename, AlignmentReaderImpl reader) throws IOException {
        if (!"1.9.5-".equals(reader.getGobyVersion())) return;

        final GZIPInputStream indexStream = new GZIPInputStream(new FileInputStream(basename + ".index"));

        final CodedInputStream codedInput = CodedInputStream.newInstance(indexStream);
        codedInput.setSizeLimit(Integer.MAX_VALUE);
        final Alignments.AlignmentIndex index = Alignments.AlignmentIndex.parseFrom(codedInput);

        LongArrayList indexOffsets = new LongArrayList();
        LongArrayList upgradedOffsets = new LongArrayList();
        LongArrayList indexAbsolutePositions = new LongArrayList();
        LongArrayList upgradedIndexAbsolutePositions = new LongArrayList();


        for (final long offset : index.getOffsetsList()) {
            indexOffsets.add(offset);
        }
        for (final long absolutePosition : index.getAbsolutePositionsList()) {
            indexAbsolutePositions.add(absolutePosition);
        }
        // trimming is essential for the binary search to work reliably with the result of elements():
        indexAbsolutePositions.trim();
        indexOffsets.trim();
        long[] targetPositionOffsets;
        int[] targetLengths = reader.getTargetLength();

        //prepare according to new indexing scheme:
        targetPositionOffsets = new long[targetLengths.length];
        targetPositionOffsets[0] = 0;
        for (int targetIndex = 1; targetIndex < targetLengths.length; targetIndex++) {
            targetPositionOffsets[targetIndex] =
                    targetLengths[targetIndex - 1] +
                            targetPositionOffsets[targetIndex - 1];
        }
        long previousAbsolutePosition = -1;
        ProgressLogger progress = new ProgressLogger(LOG);
        progress.expectedUpdates = indexOffsets.size();
        progress.priority = Level.INFO;
        progress.start();
        // push the very first entry to the index at offset zero. This is necessary because 1.9.5 did not include
        // offset information and absolute position for the very first entry of the alignment.
        Alignments.AlignmentEntry entry = reader.next();
        previousAbsolutePosition = pushEntryToIndex(upgradedOffsets, upgradedIndexAbsolutePositions,
                targetPositionOffsets, previousAbsolutePosition, 0, entry);

        for (long indexOffset : indexOffsets) {
            // for each offset in the entries file, obtain the first entry then recode absolute position:

            entry = fetchFirstEntry(reader, indexOffset);
            previousAbsolutePosition = pushEntryToIndex(upgradedOffsets, upgradedIndexAbsolutePositions,
                    targetPositionOffsets, previousAbsolutePosition, indexOffset, entry);

            progress.lightUpdate();
        }
        progress.stop();
        writeIndex(basename, upgradedOffsets, upgradedIndexAbsolutePositions);
        upgradeHeaderVersion(basename);
        if (verbose) {
            System.out.printf("alignment %s upgraded successfully.%n", basename);
        }
    }

    private long pushEntryToIndex(LongArrayList upgradedOffsets, LongArrayList upgradedIndexAbsolutePositions, long[] targetPositionOffsets, long previousAbsolutePosition, long indexOffset, Alignments.AlignmentEntry entry) {
        //   System.out.printf("entry target=%d position=%d %n", entry.getTargetIndex(), entry.getPosition());
        if (entry == null) {

            if (verbose) {
                System.err.println("Error: Cannot obtain entry at start of chunk for indexOffset: " + indexOffset);
                System.exit(10);
            }
        } else {
            int targetIndex = entry.getTargetIndex();
            int position = entry.getPosition();

            long newAbsolutePosition = targetPositionOffsets[targetIndex] + position;
            if (newAbsolutePosition > previousAbsolutePosition) {
                upgradedOffsets.add(indexOffset);
                upgradedIndexAbsolutePositions.add(newAbsolutePosition);
                previousAbsolutePosition = newAbsolutePosition;

            }
        }
        return previousAbsolutePosition;
    }

    private void upgradeHeaderVersion(String basename) throws IOException {
        InputStream headerStream;
        try {
            headerStream = new GZIPInputStream(new FileInputStream(basename + ".header"));
        } catch (IOException e) {
            // try not compressed for compatibility with 1.4-:
            LOG.trace("falling back to legacy 1.4- uncompressed header.");

            headerStream = new FileInputStream(basename + ".header");
        }
        // accept very large header messages, since these may contain query identifiers:
        final CodedInputStream codedInput = CodedInputStream.newInstance(headerStream);
        codedInput.setSizeLimit(Integer.MAX_VALUE);
        final Alignments.AlignmentHeader header = Alignments.AlignmentHeader.parseFrom(codedInput);

        Alignments.AlignmentHeader.Builder upgradedHeader = Alignments.AlignmentHeader.newBuilder(header);
        upgradedHeader.setVersion(VersionUtils.getImplementationVersion(UpgradeTo1_9_6.class));
        FileUtils.moveFile(new File(basename + ".header"), new File(basename + ".header.bak"));
        GZIPOutputStream headerOutput = new GZIPOutputStream(new FileOutputStream(basename + ".header"));
        try {
            upgradedHeader.build().writeTo(headerOutput);
        } finally {
            headerOutput.close();
        }
    }

    private void writeIndex(String basename, LongArrayList indexOffsets, LongArrayList indexAbsolutePositions) throws IOException {
        GZIPOutputStream indexOutput = null;
        try {
            FileUtils.moveFile(new File(basename + ".index"), new File(basename + ".index.bak"));
            indexOutput = new GZIPOutputStream(new FileOutputStream(basename + ".index"));
            final Alignments.AlignmentIndex.Builder indexBuilder = Alignments.AlignmentIndex.newBuilder();
            assert (indexOffsets.size() == indexAbsolutePositions.size()) : "index sizes must be consistent.";
            indexBuilder.addAllOffsets(indexOffsets);
            indexBuilder.addAllAbsolutePositions(indexAbsolutePositions);
            indexBuilder.build().writeTo(indexOutput);
        } finally {
            if (indexOutput != null) indexOutput.close();

        }
    }



    private Alignments.AlignmentEntry fetchFirstEntry(AlignmentReaderImpl reader, long indexOffset) throws IOException {
        reader.seek(indexOffset);
        if (reader.hasNext()) return reader.next();
        else return null;
    }

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(AlignmentReaderImpl.class);

    public void setSilent(boolean silent) {
        this.verbose = !silent;
    }
}
