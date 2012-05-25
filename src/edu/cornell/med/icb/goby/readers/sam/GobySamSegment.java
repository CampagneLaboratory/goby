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

import com.google.protobuf.ByteString;
import edu.cornell.med.icb.goby.util.pool.Resettable;
import edu.cornell.med.icb.goby.util.pool.ResettableObjectPoolInterface;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.bytes.ByteList;
import it.unimi.dsi.lang.MutableString;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

/**
 * Class to assist with parsing SAM to Goby, making SAM segments, one per splice piece.
 */
public class GobySamSegment implements Resettable {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(GobySamSegment.class);

    private ResettableObjectPoolInterface<GobyQuickSeqvar> gobyQuickSeqvarPool;

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
    int queryAlignedLength;
    int targetAlignedLength;
    boolean reverseStrand;

    private MutableString debugMessage;
    ByteArrayList softClippedQualityLeft;
    ByteArrayList softClippedQualityRight;

    /**
     * Not for use.
     */
    private GobySamSegment() {
    }

    public GobySamSegment(final ResettableObjectPoolInterface<GobyQuickSeqvar> gobyQuickSeqvarPool) {
        this();
        this.gobyQuickSeqvarPool = gobyQuickSeqvarPool;
        softClippedBasesLeft = new MutableString();
        softClippedBasesRight = new MutableString();
        softClippedQualityLeft = new ByteArrayList();
        softClippedQualityRight = new ByteArrayList();
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
        firstReadIndex = 1;
        firstRefPosition = 0;
        queryPosition = 0;

        distanceToNextSegment = 0;
        position = 0;
        softClippedBasesLeft.length(0);
        softClippedBasesRight.length(0);
        softClippedQualityRight.size(0);
        softClippedQualityLeft.size(0);
        for (final GobyQuickSeqvar seqvar : sequenceVariations) {
            gobyQuickSeqvarPool.returnObject(seqvar);
        }
        sequenceVariations.clear();
        readBases.length(0);
        readQuals.clear();
        refBases.length(0);
        diffBases.length(0);
        queryAlignedLength = 0;
        targetAlignedLength = 0;
        reverseStrand = false;
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
            debugMessage.append("queryPos    =").append(queryPosition).append('\n');
            debugMessage.append("qAlignLen   =").append(queryAlignedLength).append('\n');
            debugMessage.append("tAlignLen   =").append(targetAlignedLength).append('\n');
            debugMessage.append("refPosition =").append(firstRefPosition).append('\n');
            debugMessage.append("readIndex   =").append(firstReadIndex).append('\n');
            debugMessage.append("refBases    =").append(refBases.toString()).append('\n');
            debugMessage.append("readBases   =").append(readBases.toString()).append('\n');
            debugMessage.append("readQuals   =");
            debugOutputQuals().append('\n');

            for (int i = 0; i < sizeRefBases; i++) {
                final char refChar = refBases.charAt(i);
                final char readChar = readBases.charAt(i);
                if (refChar == readChar) {
                    diffBases.append('_');
                } else {
                    diffBases.append('X');
                }
            }
            debugMessage.append("diff        =").append(diffBases.toString()).append('\n');
            for (final GobyQuickSeqvar seqvar : sequenceVariations) {
                debugMessage.append(seqvar);
            }
        }
        LOG.debug(debugMessage);
    }

    public MutableString debugOutputQuals() {
        int i = 0;
        for (final Byte readQual : readQuals) {
            if (i > 0) {
                debugMessage.append(',');
            }
            debugMessage.append('[').append(i++).append(':').append(readQual).append(']');
        }
        return debugMessage;
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

    public MutableString getReadBases() {
        return readBases;
    }

    public ByteList getReadQuals() {
        return readQuals;
    }

    public MutableString getRefBases() {
        return refBases;
    }

    public boolean isReverseStrand() {
        return reverseStrand;
    }

    public int getSequenceVariationsCount() {
        return sequenceVariations.size();
    }

    public List<GobyQuickSeqvar> getSequenceVariations() {
        return sequenceVariations;
    }

    public GobyQuickSeqvar getSequenceVariations(final int i) {
        return sequenceVariations.get(i);
    }

    public void observeVariations() {
        final boolean hasReadQuals = !readQuals.isEmpty();
        GobyQuickSeqvar seqvar = null;
        int currentReadIndex = firstReadIndex - (reverseStrand ? -1 : 1);
        int currentRefPosition = firstRefPosition - 1;
        for (int i = 0; i < readBases.length(); i++) {
            final char refChar = refBases.charAt(i);
            final char readChar = readBases.charAt(i);
            if (readChar != '-') {
                // Deletion
                currentReadIndex += (reverseStrand ? -1 : 1) * 1;
            }
            if (refChar != '-') {
                // Insertion
                currentRefPosition += 1;
            }
            if (readChar == refChar) {
                continue;
            }

            boolean makeNewSeqvar = false;
            if (seqvar == null || seqvar.lastIndexPosition != i - 1) {
                makeNewSeqvar = true;
            } else {
                final int toLength = seqvar.to.length();
                final char prevReadChar = seqvar.to.charAt(toLength - 1);
                final char prevRefChar = seqvar.from.charAt(toLength - 1);
                // Don't group mutations and indels
                if ((prevReadChar == '-' && readChar != '-') || (readChar == '-' && prevReadChar != '-') ||
                     (prevRefChar == '-' && refChar != '-') || (refChar == '-' && prevRefChar != '-')) {
                    makeNewSeqvar = true;
                }
            }

            if (makeNewSeqvar) {
                seqvar = gobyQuickSeqvarPool.borrowObject();
                seqvar.lastIndexPosition = i;
                seqvar.setReadIndex(readChar == '-' && reverseStrand ? currentReadIndex - 1 : currentReadIndex);
                seqvar.setPosition(currentRefPosition);
                sequenceVariations.add(seqvar);
            } else {
                seqvar.lastIndexPosition = i;
            }
            seqvar.from.append(refChar);
            seqvar.to.append(readChar);
            if (hasReadQuals && readChar != '-') {
                seqvar.toQuals.add(readQuals.get(i));
            }
        }
    }

    public ByteString getSoftClippedQualityRight() {
        return ByteString.copyFrom(softClippedQualityRight.toByteArray());
    }
     public ByteString getSoftClippedQualityLeft() {
        return ByteString.copyFrom(softClippedQualityLeft.toByteArray());
    }
}
