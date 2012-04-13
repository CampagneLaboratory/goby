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

import edu.cornell.med.icb.goby.modes.SamHelper;
import edu.cornell.med.icb.goby.reads.QualityEncoding;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import it.unimi.dsi.Util;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.bytes.ByteList;
import it.unimi.dsi.fastutil.chars.CharArrayList;
import it.unimi.dsi.fastutil.chars.CharList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;

import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;

/**
 * Assist with converting Compact Alignments to other formats such as SAM.
 * TODO: Splices that span chromosome boundaries won't be emitted, at least for now.
 */
public class ExportableAlignmentEntryData {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(ExportableAlignmentEntryData.class);

    // This was the maximum acceptable value
    private static final byte UNKNOWN_MAPPING_VALUE = 93;
    private static final int SAM_REVERSE_STRAND_FLAG = 16;
    private boolean debug;

    private RandomAccessSequenceInterface genome;

    private CharList refBases;
    private CharList readBases;
    private MutableString readBasesOriginal;
    private ByteList qualities;
    private MutableString cigarString;
    private MutableString mismatchString;
    private MutableString invalidMessage;
    private boolean hasQualities;

    // These actual* values are generally empty except during testing when they can be fed in for debugging
    private CharList actualReads;
    private ByteList actualQualities;

    private boolean invalid;
    private int startClip;
    private int endClip;
    private int queryLength;
    private int queryAlignedLength;
    private int targetAlignedLength;
    private boolean reverseStrand;
    private QualityEncoding qualityEncoding;
    private int endTargetPositionZeroBased;

    private Alignments.AlignmentEntry alignmentEntry;
    private ReadOriginInfo readOriginInfo;
    private boolean hasReadGroups;
    private String readGroup;

    LinkedList<Integer> deleteQualityIndexes;

    /**
     * Marked private so it won't be used, always needs the genome version.
     */
    private ExportableAlignmentEntryData() {
    }

    /**
     * Constructor
     *
     * @param genome the genome accessor.
     */
    public ExportableAlignmentEntryData(final RandomAccessSequenceInterface genome,
                                        final QualityEncoding qualityEncoding) {
        debug = Util.log4JIsConfigured();
        this.genome = genome;
        this.qualityEncoding = qualityEncoding;

        refBases = new CharArrayList();
        readBases = new CharArrayList();
        readBasesOriginal = new MutableString();
        qualities = new ByteArrayList();
        cigarString = new MutableString();
        mismatchString = new MutableString();
        invalidMessage = new MutableString();

        actualReads = new CharArrayList();
        actualQualities = new ByteArrayList();

        deleteQualityIndexes = new LinkedList<Integer>();

        reset();
    }

    /**
     * Reset the state of objects as they same objects are used for every iteration so we don't have to
     * create and destroy tons of objects.
     */
    private void reset() {
        refBases.clear();
        readBases.clear();
        readBasesOriginal.setLength(0);
        qualities.clear();
        cigarString.setLength(0);
        mismatchString.setLength(0);
        invalidMessage.setLength(0);

        actualReads.clear();
        actualQualities.clear();

        invalid = false;
        startClip = 0;
        endClip = 0;
        queryLength = 0;
        queryAlignedLength = 0;
        targetAlignedLength = 0;
        reverseStrand = true;
        hasQualities = false;
        endTargetPositionZeroBased = 0;

        alignmentEntry = null;

        deleteQualityIndexes.clear();
    }

    /**
     * For spliced alignments, it is necessary to output all of the fragments at once. This will
     * duplicate this object so it can be saved for future output.
     *
     * @param from the ExportableAlignmentEntryData to duplicate
     * @return the duplicated object
     */
    public static ExportableAlignmentEntryData duplicateFrom(final ExportableAlignmentEntryData from) {
        final ExportableAlignmentEntryData to = new ExportableAlignmentEntryData(from.genome, from.qualityEncoding);

        to.refBases.addAll(from.refBases);
        to.readBases.addAll(from.readBases);
        to.readBasesOriginal.append(from.readBasesOriginal);
        to.qualities.addAll(from.qualities);
        to.cigarString.append(from.cigarString);
        to.mismatchString.append(from.mismatchString);
        to.invalidMessage.append(from.invalidMessage);

        to.actualReads.addAll(from.actualReads);
        to.actualQualities.addAll(from.actualQualities);

        to.invalid = from.invalid;
        to.startClip = from.startClip;
        to.endClip = from.endClip;
        to.queryLength = from.queryLength;
        to.queryAlignedLength = from.queryAlignedLength;
        to.targetAlignedLength = from.targetAlignedLength;
        to.reverseStrand = from.reverseStrand;
        to.hasQualities = from.hasQualities;
        to.endTargetPositionZeroBased = from.endTargetPositionZeroBased;

        to.alignmentEntry = Alignments.AlignmentEntry.newBuilder(from.alignmentEntry).build();

        // No need to propagate deleteQualityIndexes

        return to;
    }

    /**
     * The target index
     *
     * @return target index
     */
    public int getTargetIndex() {
        return alignmentEntry.getTargetIndex();
    }

    public String getTargetName() {
        return genome.getReferenceName(alignmentEntry.getTargetIndex());
    }

    /**
     * The pair flags
     *
     * @return the pair flags
     */
    public int getPairFlags() {
        int flags = alignmentEntry.hasPairFlags() ? alignmentEntry.getPairFlags() : 0;
        if (alignmentEntry.getMatchingReverseStrand()) {
            flags |= SAM_REVERSE_STRAND_FLAG;
        }
        return flags;
    }

    /**
     * Get the read name. This is always the query index converted to a String, but SAM wants a string.
     *
     * @return the read name.
     */
    public String getReadName() {
        return String.valueOf(alignmentEntry.getQueryIndex());
    }

    /**
     * Get the read name. This is always the query index converted to a String, but SAM wants a string.
     *
     * @return the read name.
     */
    public int getQueryIndex() {
        return alignmentEntry.getQueryIndex();
    }

    /**
     * The 1-based start position of the alignment before any clipping (so clipping is considered part of the
     * alignment). This is the position of the first base of the actual read as aligned.
     *
     * @return the read bases, containing "-"s for insertions.
     */
    public int getStartPosition() {
        return alignmentEntry.getPosition() + 1;
    }

    /**
     * The mapping quality of the read.
     *
     * @return mapping quality of the read.
     */
    public int getMappingQuality() {
        return alignmentEntry.getMappingQuality();
    }

    /**
     * If the entry has a pair.
     *
     * @return If the entry has a pair.
     */
    public boolean hasMate() {
        return alignmentEntry.hasPairAlignmentLink();
    }

    /**
     * Assuming hasMate() is true, get the mate's reference index (target index).
     *
     * @return If the entry has a pair.
     */
    public int getMateReferenceIndex() {
        return alignmentEntry.getPairAlignmentLink().getTargetIndex();
    }

    public int getMateAlignmentStart() {
        return alignmentEntry.getPairAlignmentLink().getPosition() + 1;
    }

    public int getInferredInsertSize() {
        return alignmentEntry.getInsertSize();
    }

    public Alignments.AlignmentEntry getAlignmentEntry() {
        return alignmentEntry;
    }

    /**
     * If the conversion from AlignmentEntry was invalid.
     *
     * @return if this object is invalid
     */
    public boolean isInvalid() {
        return invalid;
    }

    /**
     * The message describing why this object is invalid.
     *
     * @return message describing why this object is invalid.
     */
    public String getInvalidMessage() {
        return invalidMessage.toString();
    }

    /**
     * Return the read bases, which include "-"s if the alignment had deletions.
     *
     * @return the read bases, containing "-"s for insertions.
     */
    public CharList getReadBases() {
        return readBases;
    }

    /**
     * Return the ORIGINAL read bases, contains no "-" even if the alignment had deletions.
     * For clipped left or right bases, this read may contains "N"s that weren't in the original
     * read but the actual bases are unobtainable without the original reads file.
     *
     * @return The original query bases (as close as possible)
     */
    public String getReadBasesOriginal() {
        return readBasesOriginal.toString();
    }


    /**
     * Return the reference bases, which include "-"'s if the alignment had inserts.
     *
     * @return The reference bases
     */
    public CharList getReferenceBases() {
        return refBases;
    }

    /**
     * Value for CIGAR (such as for SAM).
     *
     * @return the calculated cigar string.
     */
    public String getCigarString() {
        return cigarString.toString();
    }

    /**
     * Value for MD:Z for SAM.
     *
     * @return the calculated mismatch string.
     */
    public String getMismatchString() {
        return mismatchString.toString();
    }

    /**
     * Return the readQualities. These will mostly be of value UNKNOWN_MAPPING_VALUE except for values that
     * come from SequenceVariations. When the read_quality_scores field has been populated, we can regenerate
     * all the read quality scores.
     *
     * @return the read qualities.
     */
    public ByteList getReadQualities() {
        return qualities;
    }

    /**
     * Reverse complement a base. Likely only used during tests.
     *
     * @param base the base to complement
     * @return the complemented base
     */
    private char complement(final char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'C':
                return 'G';
            case 'T':
                return 'A';
            case 'G':
                return 'C';
            default:
                return '?';
        }
    }

    /**
     * Transfer actual reads into this object. This is likely only used during tests.
     *
     * @param reads         the reads in the order of the reads file
     * @param reverseStrand if this alignment entry is reverse strand
     */
    private void setActualReads(final CharList reads, final boolean reverseStrand) {
        if (reads == null || reads.isEmpty()) {
            return;
        }
        final int size = reads.size();
        if (reverseStrand) {
            for (int i = size - 1; i >= 0; i--) {
                actualReads.add(complement(reads.get(i)));
            }
        } else {
            for (final char read : reads) {
                actualReads.add(read);
            }
        }
    }

    /**
     * Transfer actual quals into this object. This is likely only used during tests.
     *
     * @param quals         the quals in the order of the reads file
     * @param reverseStrand if this alignment entry is reverse strand
     */
    private void setActualQuals(final ByteList quals, final boolean reverseStrand) {
        if (quals == null || quals.isEmpty()) {
            return;
        }
        final int size = quals.size();
        if (reverseStrand) {
            for (int i = size - 1; i >= 0; i--) {
                actualQualities.add(quals.get(i));
            }
        } else {
            for (final byte qual : quals) {
                actualQualities.add(qual);
            }
        }
    }

    /**
     * Build an ExportableAlignmentEntryData object from an alignment entry.
     *
     * @param alignmentEntry a Goby alignment entry
     */
    public void buildFrom(final Alignments.AlignmentEntry alignmentEntry) {

        buildFrom(alignmentEntry, null,
                // some alignments will store the quality scores for the original read:
                alignmentEntry.hasReadQualityScores() ? ByteArrayList.wrap(
                        alignmentEntry.getReadQualityScores().toByteArray()) : null);
    }

    /**
     * Build an ExportableAlignmentEntryData object from an alignment entry.
     *
     * @param alignmentEntry     a Goby alignment entry
     * @param actualReadsSrc     the actual reads from the original reads file. Only provided during tests.
     * @param actualQualitiesSrc the actual qualities from the original reads file. Only provided during tests.
     */
    public void buildFrom(final Alignments.AlignmentEntry alignmentEntry,
                          final CharList actualReadsSrc, final ByteList actualQualitiesSrc) {
        reset();

        reverseStrand = alignmentEntry.getMatchingReverseStrand();
        setActualReads(actualReadsSrc, reverseStrand);
        setActualQuals(actualQualitiesSrc, reverseStrand);

        startClip = alignmentEntry.getQueryPosition();
        queryLength = alignmentEntry.getQueryLength();
        queryAlignedLength = alignmentEntry.getQueryAlignedLength();
        targetAlignedLength = alignmentEntry.getTargetAlignedLength();
        endClip = queryLength - queryAlignedLength - startClip;
        final int startPosition = alignmentEntry.getPosition();

        this.alignmentEntry = alignmentEntry;

        if (hasReadGroups) {
            final int readOriginIndex = alignmentEntry.getReadOriginIndex();
            final Alignments.ReadOriginInfo info = readOriginInfo.getInfo(readOriginIndex);
            if (info == null) {
                invalid = true;
                invalidMessage.append(String.format(
                        "Cannot export read group index=%d. Index was not found in the goby header.",
                        readOriginIndex));
                if (debug) {
                    LOG.debug(invalidMessage.toString());
                }
            } else {
                readGroup=info.getOriginId();
            }
        }
        // First obtain the number of indels
        int numInserts = 0;
        for (final Alignments.SequenceVariation seqvar : alignmentEntry.getSequenceVariationsList()) {
            final String froms = seqvar.getFrom();  // references bases. '-' means INSERTION
            final String tos = seqvar.getTo();      // read bases. '-' means INSERTION
            final int fromsLength = froms.length();
            final int tosLength = tos.length();
            if (fromsLength != tosLength) {
                invalid = true;
                invalidMessage.append("Error: Sequence variation for queryIndex=").
                        append(alignmentEntry.getQueryIndex()).
                        append(" Has an invalid sequence variation. from.length != to.length");
                if (debug) {
                    LOG.debug(invalidMessage.toString());
                }
                return;
            }
            for (int i = 0; i < froms.length(); i++) {
                final char from = froms.charAt(i);
                final char to = tos.charAt(i);
                if (from == '-' && to == '-') {
                    invalid = true;
                    invalidMessage.append("Error: Sequence variation for queryIndex=").
                            append(alignmentEntry.getQueryIndex()).
                            append(" invalid. Both 'from' and 'to' bases both equal '-'");
                    if (debug) {
                        LOG.debug(invalidMessage.toString());
                    }
                    return;
                }
                if (from == '-') {
                    numInserts += 1;
                }
            }
        }

        final boolean predefinedQuals = alignmentEntry.hasReadQualityScores();

        // Construct read & ref before any sequence variations (indels, mutations)
        final int endOfLoop = targetAlignedLength + startClip + endClip + numInserts;

        final int targetIndex = alignmentEntry.getTargetIndex();
        final char[] predefStartClips = alignmentEntry.hasSoftClippedBasesLeft() ?
                        alignmentEntry.getSoftClippedBasesLeft().toCharArray() : null;
        final char[] predefEndClips = alignmentEntry.hasSoftClippedBasesRight() ?
                        alignmentEntry.getSoftClippedBasesRight().toCharArray() : null;

        for (int i = 0; i < endOfLoop; i++) {
            final int genomePosition = i + startPosition - startClip;
            final char base = genomePosition >= 0 ? genome.get(targetIndex, i + startPosition - startClip) : 'N';
            if (i < startClip) {
                // Clipped read bases. We cannot reconstruct them, oh well.
                if (predefStartClips != null) {
                    final char clipBase = predefStartClips[i];
                    if (clipBase == '=') {
                        readBases.add(base);
                    } else {
                        readBases.add(clipBase);
                    }
                } else {
                    readBases.add('N');
                }
            } else {
                readBases.add(base);
            }
            refBases.add(base);

            // By default, we don't know qualities. Phred score of 127?
            if (!predefinedQuals) {
                qualities.add(UNKNOWN_MAPPING_VALUE);  // SAMRecord.UNKNOWN_MAPPING_VALUE is 255, which isn't a byte
            }
        }
        if (predefinedQuals) {
            for (final byte value : alignmentEntry.getReadQualityScores().toByteArray()) {
                qualities.add(value);
            }
            hasQualities = true;
        }


        // Process the seqvars
        if (!alignmentEntry.getSequenceVariationsList().isEmpty()) {
            if (debug) {
                LOG.debug("Before sequence variation:");
                LOG.debug("\n" + toString());
            }
        }
        for (final Alignments.SequenceVariation seqvar : alignmentEntry.getSequenceVariationsList()) {
            if (debug) {
                LOG.debug(seqVarToString(seqvar));
            }
            final String froms = seqvar.getFrom();  // references bases. '-' means INSERTION
            final String tos = seqvar.getTo();      // read bases. '-' means INSERTION
            final int startRefPosition = seqvar.getPosition();   // refPosition, 1-based, always numbered from left
            final byte[] toQuals = seqvar.hasToQuality() ? seqvar.getToQuality().toByteArray() : null;
            for (int i = 0; i < froms.length(); i++) {
                final char from = froms.charAt(i);
                final char to = tos.charAt(i);
                final Byte toQual = toQuals == null ? null :
                        (byte) qualityEncoding.phredQualityScoreToAsciiEncoding(toQuals[i]);

                final int refPosition = startRefPosition + i - 1; // Convert back to 0-based for list access
                if (from == '-') {
                    // Insertion, missing base in the reference.
                    refBases.add(refPosition + 1, from);
                    readBases.add(refPosition + 1, to);
                    if (toQual != null) {
                        qualities.add(refPosition + 1, toQual);
                        hasQualities = true;
                    }
                } else if (to == '-') {
                    // Deletion. Missing base in the read, but we
                    if (refBases.get(refPosition) != from) {
                        invalid = false;
                        invalidMessage.append("Error: (Deletion) Sequence variation for queryIndex=").
                                append(alignmentEntry.getQueryIndex()).
                                append(" invalid. 'from' base doesn't match actual reference base. From=").
                                append(from).append(" actual=").append(refBases.get(refPosition)).append('\n');
                        if (debug) {
                            LOG.debug(invalidMessage.toString());
                        }
                        if (invalid) {
                            return;
                        }
                    }
                    readBases.set(refPosition, to);
                    deleteQualityIndexes.addFirst(refPosition);
                } else {
                    // Mutation
                    if (refBases.get(refPosition) != from) {
                        invalid = false;
                        invalidMessage.append("Error: (Mutation) Sequence variation for queryIndex=").
                                append(alignmentEntry.getQueryIndex()).
                                append(" invalid. 'from' base doesn't match actual reference base. From=").
                                append(from).append(" actual=").append(refBases.get(refPosition)).append('\n');
                        if (debug) {
                            LOG.debug(invalidMessage.toString());
                        }
                        if (invalid) {
                            return;
                        }
                    }
                    readBases.set(refPosition, to);
                    if (toQual != null) {
                        qualities.set(refPosition, toQual);
                        hasQualities = true;
                    }
                }
            }
        }
        for (final int deleteQualityIndex : deleteQualityIndexes) {
            qualities.remove(deleteQualityIndex);
        }

        if (numInserts > 0) {
            // Inserts, clip bases to the right so we don't go beyond read length
            refBases.size(refBases.size() - numInserts);
            readBases.size(readBases.size() - numInserts);
            qualities.size(qualities.size() - numInserts);
        }
        if (endClip > 0) {
            // endClip, mark endClip number of bases to the right as N, we don't know their actual value
            final int readSize = readBases.size();
            for (int i = 0; i < endClip; i++) {
                char clipBase;
                if (predefEndClips == null) {
                    readBases.set(readSize - i - 1, 'N');
                } else {
                    readBases.set(readSize - i - 1, predefEndClips[i]);
                }
            }
        }
        for (final char readBase : readBases) {
            if (readBase != '-') {
                readBasesOriginal.append(readBase);
            }
        }
        observeReadRefDifferences();
        endTargetPositionZeroBased = alignmentEntry.getPosition() + startClip + targetAlignedLength;
        if (debug) {
            LOG.debug("\n" + toString());
        }
    }

    public List<String> getBamAttributesList() {
        return alignmentEntry.getBamAttributesList();
    }

    public void setReadGroupInfo(final ReadOriginInfo readOriginInfo) {
        this.readOriginInfo = readOriginInfo;
        hasReadGroups = true;
    }

    public String getReadGroup() {
        return readGroup;
    }

    private enum MismatchType {
        MATCH,
        INSERTION,
        DELETION,
        MISMATCH,
    }

    private enum CigarType {
        CLIP('S'),
        MATCH('M'),
        INSERTION('I'),
        DELETION('D');
        final char code;

        private CigarType(final char code) {
            this.code = code;
        }
    }

    /**
     * Merge the splice fragments into a single ExportableAlignmentEntryData. Fragments should already be in
     * their correct order based on the underlying goby alignment's next/previous splice fragments.
     * @param fragments splice fragments to merge.
     * @return the merge of fragments
     */
    public static ExportableAlignmentEntryData mergeSpliceFragments(
            final List<ExportableAlignmentEntryData> fragments) {
        final ExportableAlignmentEntryData merged = fragments.get(0);
        final Alignments.AlignmentEntry mergedAlignment = merged.alignmentEntry;
        final int size = fragments.size();
        for (int i = 1; i < size; i++) {
            final ExportableAlignmentEntryData fragment = fragments.get(i);
            final Alignments.AlignmentEntry fragmentAlignment = fragment.alignmentEntry;
            if (mergedAlignment.getTargetIndex() != fragmentAlignment.getTargetIndex()) {
                // Not outputting splices across chromosome boundaries, at least for now.
                merged.invalid = true;
                merged.invalidMessage.append("Error: Splice segments for queryIndex=").
                        append(mergedAlignment.getQueryIndex()).
                        append(" Are on different chromosomes.");
                if (merged.debug) {
                    LOG.debug(merged.invalidMessage.toString());
                }
                return merged;
            }

            merged.mergeCigars(fragment);
            if (merged.invalid) {
                // Quit early, we've detected a problem.
                return merged;
            }
            SamHelper.appendMismatches(merged.mismatchString, fragment.mismatchString);
            merged.readBasesOriginal.append(fragment.readBasesOriginal);
            merged.qualities.addAll(fragment.qualities);
        }
        return merged;
    }

    private void mergeCigars(final ExportableAlignmentEntryData other) {
        final int gap;
        final Alignments.AlignmentEntry otherAlignmentEntry = other.alignmentEntry;
        final boolean spliceInReverse;
        if (alignmentEntry.getPosition() < otherAlignmentEntry.getPosition()) {
            // this earlier than other, this cigar ends in S, other cigar starts in S
            spliceInReverse = false;
            gap = otherAlignmentEntry.getPosition() - endTargetPositionZeroBased;
        } else {
            spliceInReverse = true;
            gap = alignmentEntry.getPosition() - other.endTargetPositionZeroBased;
            swapFieldsForReverseSplice(other);
        }
        final GobyCigarElement thisLastCigar = removeLastCigar();
        final GobyCigarElement otherFirstCigar = other.removeFirstCigar();
        if (thisLastCigar.code == 'S' && otherFirstCigar.code == 'S') {
            // As expected
            readBasesOriginal.length(queryLength - thisLastCigar.size);
            other.readBasesOriginal.delete(0, otherFirstCigar.size);
            qualities.size(queryLength - thisLastCigar.size);
            other.qualities.removeElements(0, otherFirstCigar.size);
            cigarString.append(gap).append('N');
            cigarString.append(other.cigarString);
        } else {
            invalid = true;
            invalidMessage.append("Splice cigar codes were incorrect for qi=").
                    append(alignmentEntry.getQueryIndex());
        }
    }

    private void swapFieldsForReverseSplice(final ExportableAlignmentEntryData other) {
        final MutableString tempCigarString = cigarString;
        cigarString = other.cigarString;
        other.cigarString = tempCigarString;

        final MutableString tempMismatchString = mismatchString;
        mismatchString = other.mismatchString;
        other.mismatchString = tempMismatchString;

        final MutableString tempReadBasesOriginal = readBasesOriginal;
        readBasesOriginal = other.readBasesOriginal;
        other.readBasesOriginal = tempReadBasesOriginal;

        final ByteList tempQualities = qualities;
        qualities = other.qualities;
        other.qualities = tempQualities;
    }

    public GobyCigarElement removeFirstCigar() {
        final Matcher matcher = SamHelper.FIRST_CIGAR_PATTERN.matcher(cigarString);
        GobyCigarElement cigarElement = null;
        if (matcher.find()) {
            final String matchStr = matcher.group(1);
            final char matchCode = matcher.group(2).charAt(0);
            cigarElement = new GobyCigarElement(matchStr, matchCode);
            cigarString.delete(0, cigarElement.cigarWidth);
        }
        return cigarElement;
    }

    /**
     * Obtain the last cigar element.
     * @return a CigarElement
     */
    public GobyCigarElement removeLastCigar() {
        final Matcher matcher = SamHelper.LAST_CIGAR_PATTERN.matcher(cigarString);
        GobyCigarElement cigarElement = null;
        if (matcher.find()) {
            final String matchStr = matcher.group(1);
            final char matchCode = matcher.group(2).charAt(0);
            cigarElement = new GobyCigarElement(matchStr, matchCode);
            cigarString.length(cigarString.length() - cigarElement.cigarWidth);
        }
        return cigarElement;
    }

    /**
     * Now that we have constructed the ref and read bases including seqvars, construct the cigarString
     * and misamtchString (MD:Z for SAM).
     * for this alignment.
     */
    private void observeReadRefDifferences() {
        if (startClip > 0) {
            cigarString.append(startClip).append(CigarType.CLIP.code);
        }
        final int startPos = startClip;
        final int endPos = readBases.size() - endClip;
        CigarType lastCigarType = CigarType.MATCH;
        int numLastCigarType = 0;
        MismatchType lastMismatchType = MismatchType.MATCH;
        int numLastMismatchType = 0;
        for (int i = startPos; i < endPos; i++) {
            final char readBase = readBases.get(i);
            final char refBase = refBases.get(i);
            final CigarType curCigarType;
            final MismatchType curMismatchType;
            if (readBase == '-') {
                curCigarType = CigarType.DELETION;
                curMismatchType = MismatchType.DELETION;
            } else if (refBase == '-') {
                curCigarType = CigarType.INSERTION;
                curMismatchType = MismatchType.MATCH;
            } else {
                curCigarType = CigarType.MATCH;
                if (readBase == refBase) {
                    curMismatchType = MismatchType.MATCH;
                } else {
                    curMismatchType = MismatchType.MISMATCH;
                }
            }
            if (curCigarType == lastCigarType) {
                numLastCigarType++;
            } else {
                if (numLastCigarType > 0) {
                    cigarString.append(numLastCigarType).append(lastCigarType.code);
                }
                numLastCigarType = 1;
                lastCigarType = curCigarType;
            }

            if (curMismatchType != lastMismatchType) {
                if (lastMismatchType == MismatchType.MATCH && numLastMismatchType > 0) {
                    mismatchString.append(numLastMismatchType);
                }
                numLastMismatchType = 1;
                lastMismatchType = curMismatchType;
            } else {
                if (curCigarType != CigarType.INSERTION) {
                    // Mismatch=MATCH ++ CigarType=INSERTION is like a match, but we don't advance numLastMismatch
                    numLastMismatchType++;
                }
            }

            if (curMismatchType == MismatchType.DELETION && numLastMismatchType == 1) {
                mismatchString.append("^");
            }
            if (curMismatchType == MismatchType.MISMATCH || curMismatchType == MismatchType.DELETION) {
                mismatchString.append(refBase);
            }
        }
        if (numLastCigarType > 0) {
            cigarString.append(numLastCigarType).append(lastCigarType.code);
        }
        if (lastMismatchType == MismatchType.MATCH && numLastMismatchType > 0) {
            mismatchString.append(numLastMismatchType);
        }
        if (endClip > 0) {
            cigarString.append(endClip).append(CigarType.CLIP.code);
        }
    }

    /**
     * For used during debugging, logging. Output the contents of a sequence variation.
     *
     * @param seqvar a sequence variation
     * @return return a description of a sequence variation as a string
     */
    public String seqVarToString(final Alignments.SequenceVariation seqvar) {
        final byte[] toQuals = seqvar.hasToQuality() ? seqvar.getToQuality().toByteArray() : null;
        return String.format("seqvar=[from/ref:%s, to/read:%s%s, readIndex:%d, position:%d]%n",
                seqvar.getFrom(), seqvar.getTo(),
                toQuals == null ? "" : " " + ArrayUtils.toString(toQuals),
                seqvar.getReadIndex(), seqvar.getPosition());
    }

    /**
     * Describe this object as a String. For use during debugging, primarily.
     *
     * @return this object as a string
     */
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append("queryIndex =").append(getQueryIndex()).append("\n");
        if (invalid) {
            sb.append("invalidMessage").append(invalidMessage.toString()).append("\n");
        }
        sb.append("targetIndex=").append(getTargetIndex()).append("\n");
        sb.append("targetName =").append(getTargetName()).append('\n');
        sb.append("fragIndex  =").append(alignmentEntry.getFragmentIndex()).append("\n");
        sb.append("startPos   =").append(getStartPosition()).append("\n");
        sb.append("startClip  =").append(startClip).append("\n");
        sb.append("endClip    =").append(endClip).append("\n");
        sb.append("queryLength=").append(queryLength).append("\n");
        sb.append("queryAlignedLength =").append(queryAlignedLength).append("\n");
        sb.append("targetAlignedLength=").append(targetAlignedLength).append("\n");
        sb.append("reverseStrand=").append(reverseStrand).append("\n");
        basesOutput(sb,
                "refBases    =", refBases,
                "readBases   =", readBases,
                "readBases.O =", readBasesOriginal,
                "actReadBases=", actualReads,
                "diff        =");
        sb.append("actQuals    =").append(qualsToStr(actualQualities)).append("\n");
        sb.append("quals       =").append(qualsToStr(qualities)).append("\n");
        sb.append("cigar       =").append(cigarString.toString()).append("\n");
        sb.append("md:z        =").append(mismatchString.toString()).append("\n");
        return sb.toString();
    }


    /**
     * Output a string of a given list of quality values.
     *
     * @param quals the list of quality values
     * @return string of the list of quality values.
     */
    private String qualsToStr(final ByteList quals) {
        if (quals == null || quals.isEmpty()) {
            return "[] size=0";
        }
        final StringBuilder sb = new StringBuilder();
        sb.append("[");
        final int qualsSize = quals.size();
        for (int i = 0; i < qualsSize; i++) {
            if (i > 0) {
                sb.append(", ");
            }
            sb.append(String.format("%03d", quals.get(i)));
        }
        sb.append("] size=").append(quals.size());
        return sb.toString();
    }

    /**
     * Output the various kinds of bases for this object.
     *
     * @param sb                      the stringbuilder to output to
     * @param prefixRefBases          prefix string for reference bases
     * @param refBases                the reference bases
     * @param prefixReadBases         prefix string for read bases
     * @param readBases               the read bases
     * @param prefixReadBasesOriginal prefix string for readBasesOriginal
     * @param readBasesOriginal       the readBasesOriginal
     * @param prefixActualReadBases   prefix for the actualReadBases (if they exist, during testing)
     * @param actualReadBases         the actualReadBases (if they exist, during testing)
     * @param diffPrefix              the prefix for marking the differences between readBases and refBases
     */
    private void basesOutput(final StringBuilder sb,
                             final String prefixRefBases, final CharList refBases,
                             final String prefixReadBases, final CharList readBases,
                             final String prefixReadBasesOriginal, final MutableString readBasesOriginal,
                             final String prefixActualReadBases, final CharList actualReadBases,
                             final String diffPrefix) {
        final int refBasesSize = refBases.size();
        final int readBasesSize = readBases.size();
        final StringBuilder diff = new StringBuilder();
        for (int i = 0; i < Math.max(refBasesSize, readBasesSize); i++) {
            final char b1 = i < refBasesSize ? refBases.get(i) : 'x';
            final char b2 = i < readBasesSize ? readBases.get(i) : 'x';
            if (b1 == b2) {
                diff.append('+');
            } else {
                diff.append('x');
            }
        }

        sb.append(prefixRefBases).append('[');
        for (final char base : refBases) {
            sb.append(base);
        }
        sb.append("] size=").append(refBases.size()).append('\n');

        sb.append(prefixReadBases).append('[');
        for (final char base : readBases) {
            sb.append(base);
        }
        sb.append("] size=").append(readBases.size()).append('\n');

        sb.append(prefixReadBasesOriginal).append('[');
        sb.append(readBasesOriginal);
        sb.append("] size=").append(readBasesOriginal.length()).append('\n');

        if (prefixActualReadBases != null) {
            sb.append(prefixActualReadBases).append('[');
            for (final char base : actualReadBases) {
                sb.append(base);
            }
            sb.append("] size=").append(actualReadBases.size()).append('\n');
        }

        if (diff.length() > 0) {
            sb.append(diffPrefix).append('[');
            sb.append(diff.toString());
            sb.append("] size=").append(diff.length()).append('\n');
        }
    }

    class GobyCigarElement {
        int size;
        char code;
        int cigarWidth;
        protected GobyCigarElement(final String sizeStr, final char code) {
            size = Integer.parseInt(sizeStr);
            this.code = code;
            cigarWidth = sizeStr.length() + 1;
        }
        public String toString() {
            return String.format("%d%c", size, code);
        }
    }
}
