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

import edu.cornell.med.icb.goby.reads.QualityEncoding;
import edu.cornell.med.icb.goby.util.pool.QueueObjectPool;
import it.unimi.dsi.Util;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.bytes.ByteList;
import it.unimi.dsi.lang.MutableString;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.SequenceUtil;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

/**
 * Class to parse a SAMRecord. This replaces SamHelper and SplicedSamHelper.
 */
public class SamRecordParser {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(SamRecordParser.class);

    private static final List<GobySamRecordEntry> EMPTY_LIST;
    static {
        EMPTY_LIST = new ArrayList<GobySamRecordEntry>(0);
    }

    private final List<GobySamRecordEntry> segments;
    private final MutableString diffBases;
    private final QueueObjectPool<GobySamRecordEntry> gobySamRecordPool;
    private final QueueObjectPool<GobyQuickSeqvar> gobyQuickSeqvarPool;
    private final ByteList allReadQuals;
    private boolean hasReadQuals;
    private QualityEncoding qualityEncoding;

    private int numRecordsProcessed;
    private int numRecordsSkipped;

    private final boolean debug;
    private final MutableString debugMessage;

    public SamRecordParser() {
        debug = Util.log4JIsConfigured();
        gobyQuickSeqvarPool = new QueueObjectPool<GobyQuickSeqvar>() {
            @Override
            public GobyQuickSeqvar makeObject() {
                return new GobyQuickSeqvar();
            }
        };
        gobySamRecordPool = new QueueObjectPool<GobySamRecordEntry>() {
            @Override
            public GobySamRecordEntry makeObject() {
                return new GobySamRecordEntry(gobyQuickSeqvarPool);
            }
        };
        qualityEncoding = QualityEncoding.SANGER;
        segments = new ArrayList<GobySamRecordEntry>(3);
        diffBases =  new MutableString();
        numRecordsProcessed = 0;
        numRecordsSkipped = 0;
        allReadQuals = new ByteArrayList();
        debugMessage = new MutableString(); // No need to reset. this will be reset for each use.
    }

    private void reset() {
        for (final GobySamRecordEntry segment : segments) {
            gobySamRecordPool.returnObject(segment);
        }
        segments.clear();
        diffBases.length(0);
        allReadQuals.clear();
        hasReadQuals = false;
    }

    public void clear() {
        gobyQuickSeqvarPool.clear();
        gobyQuickSeqvarPool.clear();
        numRecordsProcessed = 0;
        numRecordsSkipped = 0;
    }

    public int getNumRecordsSkipped() {
        return numRecordsSkipped;
    }

    public int getNumRecordsProcessed() {
        return numRecordsProcessed;
    }

    public QualityEncoding getQualityEncoding() {
        return qualityEncoding;
    }

    public void setQualityEncoding(final QualityEncoding qualityEncoding) {
        this.qualityEncoding = qualityEncoding;
    }

    /**
     * Parse a SAMRecord. If not spliced, this will return ONE segment.
     * If the alignment entry in the samRecord was spliced, this will return one GobySamRecordEntry for
     * each splice segment, if not spliced this will return a list containing a single GobySamRecordEntry.
     * @param samRecord the samRecord to parse.
     * @return list of GobySamRecordEntry ready to be converted to AlignmentEntries (and SequenceVariantions)
     */
    public List<GobySamRecordEntry> processRead(final SAMRecord samRecord) {
        if (samRecord.getReadUnmappedFlag()) {
            numRecordsSkipped++;
            return EMPTY_LIST;
        }
        reset();

        int numInserts = 0;
        int numDeletes = 0;
        final String allRefBases = new String(SequenceUtil.makeReferenceFromAlignment(samRecord, true));
        final String allReadBases = samRecord.getReadString();
        final boolean reverseStrand = samRecord.getReadNegativeStrandFlag();
        int refStringPosition = 0;
        int readStringPosition = 0;
        int newRefStringPosition = 0;
        int newReadStringPosition = 0;
        int numSoftClippedBasesLeft = 0;
        int readIndex = reverseStrand ? samRecord.getReadLength() : 1;
        int refPosition = 1;
        int readIndexDelta;
        int refPositionDelta;
        int fragmentIndex = 0;
        final int alignmentStartPosition = samRecord.getAlignmentStart();
        final List<CigarElement> cigarElementList = samRecord.getCigar().getCigarElements();
        final int numCigarElements = cigarElementList.size();
        final String readQualsString = samRecord.getBaseQualityString();
        hasReadQuals = readQualsString != null && readQualsString.length() > 0;
        if (hasReadQuals) {
            for (int i = 0; i < readQualsString.length(); i++) {
                allReadQuals.add(qualityEncoding.asciiEncodingToPhredQualityScore(readQualsString.charAt(i)));
            }
        }

        if (debug && LOG.isDebugEnabled()) {
            debugMessage.length(0).append('\n');
            debugMessage.append("------------------------------------\n");
            debugMessage.append("Read Name=").append(samRecord.getReadName()).append('\n');
            debugMessage.append("Read length=").append(samRecord.getReadLength()).append('\n');
            debugMessage.append("Cigar=").append(samRecord.getCigarString()).append('\n');
            debugMessage.append("MD:Z=").append(samRecord.getStringAttribute("MD")).append('\n');
            debugOutputRefBases(allRefBases);
            debugOutputReadBases(allReadBases, reverseStrand);
            debugOutputReadQuals(readQualsString);
            debugMessage.append("Alignment Start Position = ").append(alignmentStartPosition).append('\n');
            debugMessage.append("Reverse strand? = ").append(reverseStrand).append('\n');
            LOG.debug(debugMessage.toString());
            LOG.debug("Cigar Ops:");
        }

        GobySamRecordEntry segment = null;
        for (int cigarElementNum = 0; cigarElementNum < numCigarElements; cigarElementNum++) {
            final CigarElement cigarElement = cigarElementList.get(cigarElementNum);
            final CigarOperator cigarOperator = cigarElement.getOperator();
            final int cigarLength = cigarElement.getLength();
            if (segment == null) {
                segment = gobySamRecordPool.borrowObject();
                segments.add(segment);

                segment.fragmentIndex = fragmentIndex++;

                segment.targetIndex = samRecord.getReferenceIndex();
                if (samRecord.getReadPairedFlag() && !samRecord.getMateUnmappedFlag()) {
                    segment.hasMate = true;
                    segment.mateTargetIndex = samRecord.getMateReferenceIndex();
                    segment.mateStartPosition = samRecord.getMateAlignmentStart() - 1;
                } else {
                    segment.hasMate = false;
                }
                segment.pairFlags = samRecord.getFlags();

                segment.position = alignmentStartPosition - 1 + refStringPosition - numDeletes - numSoftClippedBasesLeft;
                segment.reverseStrand = reverseStrand;
                segment.readNum = numRecordsProcessed;
            }

            readIndexDelta = 0;
            refPositionDelta = 0;
            if (debug && LOG.isDebugEnabled()) {
                LOG.debug("--new cigar element--");
            }
            if (cigarOperator == CigarOperator.D) {
                numDeletes += cigarLength;
            } else if (cigarOperator == CigarOperator.I) {
                numInserts += cigarLength;
            }

            if (debug && LOG.isDebugEnabled()) {
                LOG.debug(String.format("   op=%s len=%d consumesReads=%s consumesRefs=%s",
                        cigarOperator.name(), cigarLength,
                        cigarOperator.consumesReadBases() ? "Yes" : "No",
                        cigarOperator.consumesReferenceBases() ? "Yes" : "No"));
            }

            if (cigarOperator.consumesReadBases()) {
                if (cigarOperator == CigarOperator.S) {
                    // Q: This will save soft clipping with leftmost cigar element and right most, is that enought?
                    if (cigarElementNum == 0) {
                        segment.softClippedBasesLeft.append(allReadBases.substring(readStringPosition, readStringPosition + cigarLength));
                        numSoftClippedBasesLeft += cigarLength;
                    } else if (cigarElementNum == numCigarElements - 1) {
                        segment.softClippedBasesRight.append(allReadBases.substring(readStringPosition, readStringPosition + cigarLength));
                    }
                } else {
                    segment.setPositions(readStringPosition, readIndex, refPosition);
                    segment.readBases.append(allReadBases.substring(readStringPosition, readStringPosition + cigarLength));
                    if (hasReadQuals) {
                        segment.readQuals.addAll(allReadQuals.subList(readStringPosition, readStringPosition + cigarLength));
                    }
                    segment.queryAlignedLength += cigarLength;
                }
                newReadStringPosition = readStringPosition + cigarLength;
                readIndexDelta += cigarLength;
            } else {
                if (cigarOperator == CigarOperator.N) {
                    segment.distanceToNextSegment = cigarLength;
                    segment = null;
                    if (debug && LOG.isDebugEnabled()) {
                        LOG.debug("Splice length = " + cigarLength);
                    }
                    // Reset refPosition back to 1 after a splice
                    refPositionDelta = -refPosition - cigarLength + 1;
                } else if (cigarOperator == CigarOperator.D) {
                    for (int i = 0; i < cigarLength; i++) {
                        segment.readBases.append("-");
                        if (hasReadQuals) {
                            segment.readQuals.add((byte) -1);
                        }
                    }
                }
            }
            if (cigarOperator.consumesReferenceBases()) {
                if (cigarOperator != CigarOperator.N) {
                    segment.refBases.append(allRefBases.substring(refStringPosition, refStringPosition + cigarLength));
                    segment.targetAlignedLength += cigarLength;
                }
                refPositionDelta += cigarLength;
            } else {
                if (cigarOperator == CigarOperator.I) {
                    for (int i = 0; i < cigarLength; i++) {
                        segment.refBases.append("-");
                    }
                }
            }
            // ALWAYS increment refPosition by cigarLength because of how we build allRefBases.
            newRefStringPosition = refStringPosition + cigarLength;

            refStringPosition = newRefStringPosition;
            readStringPosition = newReadStringPosition;

            readIndex += (reverseStrand ? -1 : 1) * readIndexDelta;
            refPosition += refPositionDelta;
        }

        for (final GobySamRecordEntry aSegment : segments) {
            aSegment.observeVariations();
            if (debug && LOG.isDebugEnabled()) {
                aSegment.debugOutput();
            }
        }

        numRecordsProcessed++;
        return segments;
    }

    /**
     * For debugging. Append to the debug string to output the reference bases along with their reference positions.
     * @param refs the reference bases.
     */
    private void debugOutputRefBases(final String refs) {
        debugMessage.append("             ");
        int refPosition = 0;
        for (final char c : refs.toCharArray()) {
            if (c == 'N') {
                refPosition = 0;
                debugMessage.append("_");
            } else {
                if (c != '-' && c != '0') {
                    refPosition++;
                }
                debugMessage.append(refPosition % 10);
            }
        }
        debugMessage.append('\n');
        debugMessage.append("allRefBases =").append(refs).append('\n');
    }

    /**
     * For debugging. Append to the debug string to output the read bases along with their read positions.
     * @param reads the read bases
     * @param reverseStrand if the read was aligned in the reverse strand
     */
    private void debugOutputReadBases(final String reads, final boolean reverseStrand) {
        debugMessage.append("allReadBases=").append(reads).append('\n');
        debugMessage.append("             ");
        if (reverseStrand) {
            for (int i = reads.length() ; i >= 1; i--) {
                debugMessage.append(i % 10);
            }
        } else {
            for (int i = 0; i < reads.length(); i++) {
                debugMessage.append((i + 1) % 10);
            }
        }
        debugMessage.append('\n');
    }

    /**
     * For debugging. Append to the debug string to output the read qualities both in the original ASCII
     * and as transformed by the qualityEncoder.
     * @param readQualsString the read qualities string, exactly as it appears in the SAM file
     * @params reverseStrand if the read was aligned in the reverse strand
     */
    private void debugOutputReadQuals(final String readQualsString) {
        debugMessage.append("readQualsStr=").append(readQualsString).append('\n');
        debugMessage.append("allReadQuals=");
        if (hasReadQuals) {
            for (int i = 0; i < allReadQuals.size(); i++) {
                if (i > 0) {
                    debugMessage.append(',');
                }
                debugMessage.append('[').append(i).append(':').append(allReadQuals.get(i)).append(']');
            }
        } else {
            debugMessage.append("None");
        }
        debugMessage.append('\n');
    }
}
