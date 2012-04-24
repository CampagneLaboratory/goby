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

import edu.cornell.med.icb.goby.util.pool.QueueObjectPool;
import edu.cornell.med.icb.goby.util.pool.Resettable;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.bytes.ByteList;
import it.unimi.dsi.lang.MutableString;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

/**
 * Class to assist with parsing SAM to Goby, making SAM segments, one per splice piece.
 * SamRecordParser, when parsing a non-spliced read, will create one of these.
 * If the read has three splices, it will create three of these.
 * AlignmentEntry objects will be created from these objects.
 */
public class GobySamRecord implements Resettable {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(GobySamRecord.class);

    private QueueObjectPool<GobyQuickSeqvar> gobyQuickSeqvarPool;

    boolean firstPositionsSet;
    int firstReadIndex;
    int firstRefPosition;
    int queryPosition;

    int distanceToNextSegment;
    int position;
    MutableString softClippedBasesLeft;
    MutableString softClippedBasesRight;
    List<GobyQuickSeqvar> sequenceVariations;
    MutableString readBases;
    ByteList readQuals;
    MutableString refBases;
    MutableString diffBases;
    /** Number of non-N parts from the cigar this string is made from. */
    int parts;
    int numInserts;
    int numDeletes;
    int queryAlignedLength;
    int targetAlignedLength;
    boolean reverseStrand;
    int fragmentIndex;
    int readNum;

    int targetIndex;
    boolean hasMate;
    int pairFlags;
    int mateTargetIndex;
    int mateStartPosition;

    private MutableString debugMessage;

    /**
     * Not for use.
     */
    private GobySamRecord() {
    }

    public GobySamRecord(final QueueObjectPool<GobyQuickSeqvar> gobyQuickSeqvarPool) {
        this();
        this.gobyQuickSeqvarPool = gobyQuickSeqvarPool;
        softClippedBasesLeft = new MutableString();
        softClippedBasesRight = new MutableString();
        sequenceVariations = new ArrayList<GobyQuickSeqvar>();
        readBases = new MutableString();
        readQuals = new ByteArrayList();
        refBases = new MutableString();
        diffBases = new MutableString();
        debugMessage = new MutableString();  // No need to reset. this will be reset for each use.
        reset();
    }

    @Override
    public void reset() {
        firstPositionsSet = false;
        firstReadIndex = 0;
        firstRefPosition = 0;
        queryPosition = 0;

        distanceToNextSegment = 0;
        position = 0;
        softClippedBasesLeft.length(0);
        softClippedBasesRight.length(0);
        for (final GobyQuickSeqvar seqvar : sequenceVariations) {
            try {
                gobyQuickSeqvarPool.returnObject(seqvar);
            } catch (Exception e) {
                // ?? What to do with this?
            }
        }
        sequenceVariations.clear();
        readBases.length(0);
        readQuals.clear();
        refBases.length(0);
        diffBases.length(0);
        parts = 0;
        numInserts = 0;
        numDeletes = 0;
        queryAlignedLength = 0;
        targetAlignedLength = 0;
        reverseStrand = false;
        fragmentIndex = 0;
        readNum = 0;

        targetIndex = 0;
        hasMate = false;
        pairFlags = 0;
        mateTargetIndex = 0;
        mateStartPosition = 0;
    }

    public void setPositions(final int queryPosition, final int readIndex, final int refPosition) {
        if (!firstPositionsSet) {
            firstPositionsSet = true;
            this.queryPosition = queryPosition;
            firstReadIndex = readIndex;
            firstRefPosition = refPosition;
        }
    }

    public void debugOutput() {
        debugMessage.length(0);
        diffBases.length(0);
        final int sizeRefBases = refBases.length();
        final int sizeReadBases = readBases.length();
        if (sizeReadBases == 0 && sizeRefBases == 0) {
            debugMessage.append("No bases. Probably padding.");
        } else {
            debugMessage.append('\n');
            debugMessage.append("alignStart  =").append(position).append('\n');
            debugMessage.append("fragIndex   =").append(fragmentIndex).append('\n');
            debugMessage.append("queryPos    =").append(queryPosition).append('\n');
            debugMessage.append("qAlignLen   =").append(queryAlignedLength).append('\n');
            debugMessage.append("tAlignLen   =").append(targetAlignedLength).append('\n');
            debugMessage.append("refPosition =").append(firstRefPosition).append('\n');
            debugMessage.append("readIndex   =").append(firstReadIndex).append('\n');
            debugMessage.append("refBases    =").append(refBases.toString()).append('\n');
            debugMessage.append("readBases   =").append(readBases.toString()).append('\n');

            for (int i = 0; i < sizeRefBases; i++) {
                final char refChar = refBases.charAt(i);
                final char readChar = readBases.charAt(i);
                if (refChar == readChar) {
                    diffBases.append('_');
                } else {
                    diffBases.append('X');
                }
            }
            debugMessage.append("diff        =").append(diffBases.toString());
        }
        LOG.debug(debugMessage);
    }

    public void observeVariations() {
        final boolean hasReadQuals = !readQuals.isEmpty();
        GobyQuickSeqvar seqvar = null;
        int currentReadIndex = firstReadIndex - (reverseStrand ? -1 : 1);
        int currentRefPosition = firstRefPosition - 1;
        int readIndexDelta;
        int refPositionDelta;
        for (int i = 0; i < readBases.length(); i++) {
            readIndexDelta = 1;
            refPositionDelta = 1;
            final char refChar = refBases.charAt(i);
            final char readChar = readBases.charAt(i);
            if (readChar == '-') {
                // Deletion
                readIndexDelta = 0;
            } else if (refChar == '-') {
                // Insertion
                refPositionDelta = 0;
            }
            currentReadIndex += (reverseStrand ? -1 : 1) * readIndexDelta;
            currentRefPosition += refPositionDelta;
            if (readChar == refChar) {
                continue;
            }

            if (seqvar == null) {
                seqvar = gobyQuickSeqvarPool.borrowObject();
                seqvar.lastIndexPosition = i;
                seqvar.readIndex = currentReadIndex;
                seqvar.position = currentRefPosition;
                sequenceVariations.add(seqvar);
            } else {
                boolean makeNewSeqvar = false;
                if (seqvar.lastIndexPosition != i - 1) {
                    makeNewSeqvar = true;
                } else if (hasReadQuals && readChar == '-') {
                    // If readChar is '-' and we have quals, we only want to append to the previous
                    // seqvar if the previous to was a '-'.
                    final int toLength = seqvar.to.length();
                    if (toLength > 0 && seqvar.to.charAt(toLength - 1) != '-') {
                        makeNewSeqvar = true;
                    }
                }

                if (makeNewSeqvar) {
                    seqvar = gobyQuickSeqvarPool.borrowObject();
                    seqvar.lastIndexPosition = i;
                    seqvar.readIndex = currentReadIndex;
                    seqvar.position = currentRefPosition;
                    sequenceVariations.add(seqvar);
                } else {
                    seqvar.lastIndexPosition = i;
                }
            }
            seqvar.from.append(refChar);
            seqvar.to.append(readChar);
            if (hasReadQuals && readChar != '-') {
                seqvar.toQuals.add(readQuals.get(i));
            }
        }
    }

    public int getPosition() {
        return position;
    }

    public int getQueryPosition() {
        return queryPosition;
    }

    public String getSoftClippedBasesLeft() {
        return softClippedBasesLeft.toString();
    }

    public String getSoftClippedBasesRight() {
        return softClippedBasesRight.toString();
    }

    public int getQueryAlignedLength() {
        return queryAlignedLength;
    }

    public int getTargetAlignedLength() {
        return targetAlignedLength;
    }

    public int getFragmentIndex() {
        return fragmentIndex;
    }

    public int getSequenceVariationsCount() {
        return sequenceVariations.size();
    }

    public GobyQuickSeqvar getSequenceVariations(final int i) {
        return sequenceVariations.get(i);
    }

}
