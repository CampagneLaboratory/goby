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

package edu.cornell.med.icb.goby.alignments;

import com.google.protobuf.CodedInputStream;
import edu.cornell.med.icb.goby.GobyVersion;
import edu.cornell.med.icb.goby.modes.ConcatenateAlignmentMode;
import edu.cornell.med.icb.goby.reads.MessageChunksReader;
import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.Date;
import java.util.UUID;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Fix the index of large alignments (>2GB). Versions prior to Goby 1.9.8.2 had an issue where the index was missing
 * entries for entries beyond the 2GB mark.
 *
 * @author Fabien Campagne
 *         Date: Jan 20, 2012
 *         Time: 2:25:21 PM
 */
public class UpgradeTo1_9_8_2 {
    private boolean verbose;

    public void upgrade(final String basename, final AlignmentReaderImpl reader) throws IOException {
        if (GobyVersion.isMoreRecent(reader.getGobyVersion(), "goby_1.9.8.2")) {
            return;
        }
        String entriesFilename = basename + ".entries";
        File entriesFile = new File(entriesFilename);
        if (entriesFile.length() < Integer.MAX_VALUE) {
            // entries file smaller than 2GB, nothing to do.
            //System.out.println("Entries file smaller than 2GB, would ignore");
            return;
        }
        System.out.printf("Upgrading %s.. This will take a while.%n", basename);
        ConcatenateAlignmentMode concat = new ConcatenateAlignmentMode();

        String[] inputs = {basename};
        concat.setInputFileNames(inputs);
        String tmpBasename = makeTmpBasename(basename);
        concat.setUpgrade(false);
        concat.setOutputFilename(tmpBasename);
        concat.setIgnoreTmh(true);
        try {
            concat.execute();

            if (new File(tmpBasename + ".index").exists()) {

                FileUtils.moveFile(new File(basename + ".index"), new File(makeBackFilename(basename + ".index", ".bak")));
                FileUtils.moveFile(new File(tmpBasename + ".index"), new File(basename + ".index"));


                upgradeHeaderVersion(basename);
                if (verbose) {
                    System.out.printf("alignment %s upgraded successfully.%n", basename);
                }
            }

        else{
            if (verbose) {
                System.out.printf("alignment %s failed to upgrade to 1.9.8.2.%n", basename);
            }
        }
    }

    finally

    {
        FileUtils.deleteQuietly(new File(tmpBasename + ".index"));
        FileUtils.deleteQuietly(new File(tmpBasename + ".entries"));
        FileUtils.deleteQuietly(new File(tmpBasename + ".tmh"));
        FileUtils.deleteQuietly(new File(tmpBasename + ".header"));
        FileUtils.deleteQuietly(new File(tmpBasename + ".stats"));
    }

}

    private String makeTmpBasename(String basename) {
        return basename + "-TMP-"+ new Date().getTime() + UUID.randomUUID();
    }

    /**
     * Returns true upon success.
     *
     * @param entriesFilename
     * @param indexOffsets
     * @param absolutePositions
     * @param countVoid
     * @return
     * @throws IOException
     */
    private boolean scanEntriesPast2G(String entriesFilename,
                                      LongArrayList indexOffsets,
                                      LongArrayList absolutePositions,
                                      int countVoid) throws IOException {

        MessageChunksReader reader = new MessageChunksReader(new FileInputStream(entriesFilename));
        int count = 0;
        int addedCount = 0;
        long position = reader.position();
        while (reader.hasNext(null, 0)) {

            if (count > indexOffsets.size()) {
                indexOffsets.add(position);
                addedCount++;
            } else {
                if (count < indexOffsets.size()) {
                    final long previousPos = indexOffsets.getLong(count);
                    if (position != 0 && position != previousPos) {
                        System.out.printf("recalculated position %d must not clash with previous position %d for <2GB part of file count=%d %n",
                                position, previousPos, count);
                    }
                } else {
                    System.out.println("past end of previous offsets: " + count);
                }
            }
            count++;
            position = reader.position();

        }

        assert addedCount == countVoid : "addedCount must match void count.";
        if (indexOffsets.size() != absolutePositions.size()) {
            System.out.printf("new sizes must match indexOffsets.size=%d absolutePositions.size=%d %n",
                    indexOffsets.size(),
                    absolutePositions.size());
            return false;
        }
        return true;
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
        upgradedHeader.setVersion(VersionUtils.getImplementationVersion(UpgradeTo1_9_8_2.class));
        FileUtils.moveFile(new File(basename + ".header"), new File(makeBackFilename(basename + ".header", ".bak")));
        GZIPOutputStream headerOutput = new GZIPOutputStream(new FileOutputStream(basename + ".header"));
        try {
            upgradedHeader.build().writeTo(headerOutput);
        } finally {
            headerOutput.close();
        }
    }

    private String makeBackFilename(String filename, String ext) {
        int counter = 1;
        String extCount = ext;
        while (new File(filename + extCount).exists()) {
            counter++;
            extCount = ext + Integer.toString(counter);
        }
        return filename + extCount;
    }

    private void writeIndex(String basename, LongArrayList indexOffsets, LongArrayList indexAbsolutePositions) throws IOException {
        GZIPOutputStream indexOutput = null;
        try {
            FileUtils.moveFile(new File(basename + ".index"), new File(makeBackFilename(basename + ".index", ".bak")));
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
