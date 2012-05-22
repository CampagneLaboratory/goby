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
import edu.cornell.med.icb.goby.util.pool.QueueResettableObjectPool;
import edu.cornell.med.icb.goby.util.pool.Resettable;
import edu.cornell.med.icb.goby.util.pool.ResettableObjectPoolInterface;
import it.unimi.dsi.Util;
import it.unimi.dsi.fastutil.bytes.ByteList;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

/**
 * A sam record. Can contain one or more segments (more than one if the sam record is spliced).
 * Each segment will become a Goby AlignmentEntry.
 */
public class GobySamRecord implements Resettable {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(GobySamRecord.class);

    private final ResettableObjectPoolInterface<GobySamSegment> gobySamSegmentPool;
    private final ResettableObjectPoolInterface<GobyQuickSeqvar> gobyQuickSeqvarPool;

    List<GobyQuickSeqvar> allSequenceVariations;

    List<GobySamSegment> segments;

    // Data owned by SamRecordParser, do not reset.
    String query;

    // Data owned by SamRecordParser, do not reset
    ByteList readQuals;

    int readNum;
    int targetIndex;
    boolean reverseStrand;
    boolean hasMate;
    int pairFlags;
    int mateTargetIndex;
    int mateStartPosition;

    int numInserts;
    int numDeletes;
    int queryAlignedLength;
    int targetAlignedLength;

    private final boolean debug;

    /**
     * Not for use.
     */
    public GobySamRecord() {
        debug = Util.log4JIsConfigured();
        gobySamSegmentPool = new QueueResettableObjectPool<GobySamSegment>() {
            @Override
            public GobySamSegment makeObject() {
                return new GobySamSegment(gobyQuickSeqvarPool);
            }
        };
        gobyQuickSeqvarPool = new QueueResettableObjectPool<GobyQuickSeqvar>() {
            @Override
            public GobyQuickSeqvar makeObject() {
                return new GobyQuickSeqvar();
            }
        };
        allSequenceVariations = new ArrayList<GobyQuickSeqvar>();
        segments = new ArrayList<GobySamSegment>();
        reset();
    }

    @Override
    public void reset() {
        for (final GobySamSegment segment : segments) {
            gobySamSegmentPool.returnObject(segment);
        }
        segments.clear();
        allSequenceVariations.clear(); // These are the same ones in segments, no need to return these to pool
        readNum = 0;
        targetIndex = 0;
        reverseStrand = false;
        hasMate = false;
        pairFlags = 0;
        mateTargetIndex = 0;
        mateStartPosition = 0;

        numInserts = 0;
        numDeletes = 0;
        queryAlignedLength = 0;
        targetAlignedLength = 0;
    }

    public int getQueryLength() {
        return query.length();
    }

    public int getNumSegments() {
        return segments.size();
    }

    public List<GobySamSegment> getSegments() {
        return segments;
    }

    public GobySamSegment getSegment(final int index) {
        return segments.get(index);
    }

    public GobySamSegment addSegment() {
        final GobySamSegment segment = gobySamSegmentPool.borrowObject();
        segments.add(segment);
        return segment;
    }

    public int getReadNum() {
        return readNum;
    }

    public int getTargetIndex() {
        return targetIndex;
    }

    public boolean isReverseStrand() {
        return reverseStrand;
    }

    public boolean isHasMate() {
        return hasMate;
    }

    public int getPairFlags() {
        return pairFlags;
    }

    public int getMateTargetIndex() {
        return mateTargetIndex;
    }

    public int getMateStartPosition() {
        return mateStartPosition;
    }

    public int getNumInserts() {
        return numInserts;
    }

    public int getNumDeletes() {
        return numDeletes;
    }

    public int getTargetAlignedLength() {
        return targetAlignedLength;
    }

    public int getQueryAlignedLength() {
        return queryAlignedLength;
    }

    public int getSequenceVariationsCount() {
        return allSequenceVariations.size();
    }

    public GobyQuickSeqvar getSequenceVariations(final int i) {
        return allSequenceVariations.get(i);
    }

    public List<GobyQuickSeqvar> getSequenceVariations() {
        return allSequenceVariations;
    }

    public String getQuery() {
        return query;
    }

    public byte[] getReadQualitiesAsBytes() {
        return readQuals.toByteArray();
    }

    public void observeVariations() {
        for (final GobySamSegment segment : segments) {
            segment.observeVariations();
            if (debug && LOG.isDebugEnabled()) {
                segment.debugOutput();
            }
            allSequenceVariations.addAll(segment.sequenceVariations);
            targetAlignedLength += segment.targetAlignedLength;
            queryAlignedLength += segment.queryAlignedLength;
        }
    }


}
