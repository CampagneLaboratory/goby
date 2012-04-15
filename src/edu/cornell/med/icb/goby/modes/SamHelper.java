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

package edu.cornell.med.icb.goby.modes;

import edu.cornell.med.icb.goby.reads.QualityEncoding;
import it.unimi.dsi.Util;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.lang.MutableString;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class to assist with parsing SAM files. Translated from C_Alignments.cc.
 *
 * @author Kevin Dorff
 */
public class SamHelper {

    private static final Pattern CIGAR_REGEX = Pattern.compile("([0-9]+)([SMID])");
    private static final Pattern MD_REGEX = Pattern.compile("([0-9]+|[ACGTN]|\\^[ACGTN]+)");
    private static final Pattern NUMERIC_REGEX = Pattern.compile("^[0-9]+$");

    private static final Pattern FIRST_NUMBER_PATTERN = Pattern.compile("^(\\d+)");
    private static final Pattern LAST_NUMBER_PATTERN = Pattern.compile("(\\d+)$");

    public static final Pattern FIRST_CIGAR_PATTERN = Pattern.compile("^(\\d+)([MIDNSHP])");
    public static final Pattern LAST_CIGAR_PATTERN = Pattern.compile("(\\d+)([MIDNSHP])$");

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(SamHelper.class);

    private byte minQualValue = 0;
    private MutableString cigar = new MutableString();
    private MutableString md = new MutableString();
    private MutableString sourceQuery = new MutableString();
    private MutableString sourceQual = new MutableString();
    private MutableString query = new MutableString();
    private MutableString qual = new MutableString();
    private MutableString ref = new MutableString();
    private int alignedLength;
    private int queryAlignedLength;
    private int targetAlignedLength;
    private int numInsertions;
    private int numDeletions;
    private int numMisMatches;
    private int score;
    private int numLeftClipped;
    private int numRightClipped;
    private int position;
    private int queryIndex;
    private int queryPosition;
    private int queryLength;
    private boolean reverseStrand;
    private IntList refPositions = new IntArrayList();
    private IntList readIndexes = new IntArrayList();
    private ObjectList<SamSequenceVariation> sequenceVariations = new ObjectArrayList<SamSequenceVariation>();

    private MutableString logval = new MutableString();

    private QualityEncoding qualityEncoding = QualityEncoding.SANGER;
    private static boolean debug = true;


    public SamHelper() {
        // don't even dare go through the debugging code if log4j was not configured. The debug code
        // is way too slow to run unintentionally in production!
        debug = debug && Util.log4JIsConfigured();
    }

    public void reset() {
        cigar.setLength(0);
        md.setLength(0);
        sourceQuery.setLength(0);
        sourceQual.setLength(0);
        query.setLength(0);
        qual.setLength(0);
        ref.setLength(0);
        alignedLength = 0;
        numInsertions = 0;
        numDeletions = 0;
        numMisMatches = 0;
        score = 0;
        numLeftClipped = 0;
        numRightClipped = 0;
        position = 0;
        queryIndex = 0;
        queryPosition = 0;
        queryLength = 0;
        reverseStrand = false;
        refPositions.clear();
        readIndexes.clear();
        sequenceVariations.clear();
    }

    public void setSource(final int queryIndex, final CharSequence sourceQuery,
                          final CharSequence sourceQual,
                          final CharSequence cigar,
                          final CharSequence md,
                          final int position,
                          final boolean reverseStrand,
                          int readLength) {
        this.queryLength = readLength;
        if (debug && LOG.isDebugEnabled()) {
            LOG.debug("------ new setSource --------------------------------");
            LOG.debug("position=" + (position - 1));
            LOG.debug("queryIndex=" + queryIndex);
        }
        reset();
        this.queryIndex = queryIndex;
        this.sourceQuery.setLength(0);
        if (sourceQuery != null) {
            this.sourceQuery.append(sourceQuery);
            queryLength = sourceQuery.length();
        }
        this.sourceQual.setLength(0);
        if (sourceQual != null) {
            this.sourceQual.append(sourceQual);
        }
        this.cigar.setLength(0);
        if (cigar != null) {
            this.cigar.append(cigar);
        }
        this.md.setLength(0);
        if (md != null) {
            this.md.append(md);
            this.md.toUpperCase();
        }
        this.position = position - 1;  // SAM positions are 1-based, goby are 0-based
        this.reverseStrand = reverseStrand;
        constructRefAndQuery();
        findSequenceVariations();
        SamSequenceVariation.merge(sequenceVariations);
    }

    public void setSourceWithReference(final int queryIndex, final SAMRecord samRecord, final String sourceRef) {

        final CharSequence sourceQuery = samRecord.getReadString();
        final CharSequence sourceQual = samRecord.getBaseQualityString();
        final int position = samRecord.getAlignmentStart();
        final boolean reverseStrand = samRecord.getReadNegativeStrandFlag();

        setSourceWithReference(queryIndex, sourceRef, sourceQuery, sourceQual, position, reverseStrand);
    }

    public void setSourceWithReference(final int queryIndex, final CharSequence sourceRef,
                                       final CharSequence sourceQuery,
                                       final CharSequence sourceQual, final int position,
                                       final boolean reverseStrand) {
        if (debug && LOG.isDebugEnabled()) {
            LOG.debug("------ new setSourceWithReference --------------------------------");
            LOG.debug("position=" + (position - 1));
            LOG.debug("queryIndex=" + queryIndex);
        }
        reset();
        this.queryIndex = queryIndex;
        this.sourceQuery.setLength(0);
        if (sourceQuery != null) {
            this.sourceQuery.append(sourceQuery);
            queryLength = sourceQuery.length();
        }
        this.sourceQual.setLength(0);
        if (sourceQual != null) {
            this.sourceQual.append(sourceQual);
        }
        ref.setLength(0);
        if (sourceRef != null) {
            ref.append(sourceRef);
        }
        query.setLength(0);
        if (sourceQuery != null) {
            query.append(sourceQuery);
        }
        qual.setLength(0);
        if (sourceQual != null) {
            qual.append(sourceQual);
        }
        this.position = position - 1;  // SAM positions are 1-based, goby are 0-based
        queryPosition = 0;
        this.reverseStrand = reverseStrand;
        this.numInsertions = 0;
        this.numDeletions = 0;
        this.numLeftClipped = 0;
        this.numRightClipped = 0;
        this.alignedLength = query.length();
        this.queryAlignedLength = query.length();
        this.queryLength = query.length();
        this.targetAlignedLength = ref.length();
        int scanLength = Math.min(query.length(), ref.length());
        this.numMisMatches = 0;
        for (int i = 0; i < scanLength; i++) {
            if (query.charAt(i) != ref.charAt(i)) {
                numMisMatches++;
                ref.setCharAt(i, Character.toLowerCase(ref.charAt(i)));
            }
        }
        this.score = alignedLength - numMisMatches;
        findSequenceVariations();
        SamSequenceVariation.merge(sequenceVariations);
        if (debug && LOG.isDebugEnabled()) {
            for (SamSequenceVariation var : sequenceVariations) {
                LOG.debug("... Variation " + var.toString());
            }
        }
    }


    public List<SamSequenceVariation> getSequenceVariations() {
        return sequenceVariations;
    }

    public MutableString getSourceQuery() {
        return sourceQuery;
    }

    public MutableString getQuery() {
        return query;
    }

    public MutableString getSourceQual() {
        return sourceQual;
    }

    public byte[] getSourceQualAsBytes() {

        final int length = sourceQual.length();
        final byte[] result = new byte[length];
        for (int i = 0; i < length; i++) {
            result[i] = qualityEncoding.asciiEncodingToPhredQualityScore(sourceQual.charAt(i));
        }
        return result;
    }

    public MutableString getQual() {
        return qual;
    }

    public MutableString getRef() {
        return ref;
    }

    public MutableString getCigar() {
        return cigar;
    }

    public MutableString getMd() {
        return md;
    }

    public void setMinQualValue(final byte minQualValue) {
        this.minQualValue = minQualValue;
    }

    public void setMinQualValue(final char minQualValue) {
        this.minQualValue = (byte) minQualValue;
    }

    public byte getMinQualValue() {
        return minQualValue;
    }

    public int getAlignedLength() {
        return alignedLength;
    }

    public int getQueryAlignedLength() {
        return queryAlignedLength;
    }

    public int getQueryLength() {
        return queryLength;
    }

    public int getTargetAlignedLength() {
        return targetAlignedLength;
    }

    public int getNumInsertions() {
        return numInsertions;
    }

    public int getNumDeletions() {
        return numDeletions;
    }

    public int getNumMisMatches() {
        return numMisMatches;
    }

    public int getScore() {
        return score;
    }

    public int getNumLeftClipped() {
        return numLeftClipped;
    }

    public int getNumRightClipped() {
        return numRightClipped;
    }

    /**
     * Return zero-based position.
     *
     * @return zero-based position.
     */
    public int getPosition() {
        return position;
    }

    public int getQueryIndex() {
        return queryIndex;
    }

    public int getQueryPosition() {
        return queryPosition;
    }

    public boolean isReverseStrand() {
        return reverseStrand;
    }

    protected void constructRefAndQuery() {
        query.setLength(0);
        qual.setLength(0);
        ref.setLength(0);
        alignedLength = 0;
        queryAlignedLength = 0;
        targetAlignedLength = 0;
        numInsertions = 0;
        numDeletions = 0;
        numMisMatches = 0;
        numLeftClipped = 0;
        numRightClipped = 0;
        score = 0;

        if (debug && LOG.isDebugEnabled()) {
            LOG.debug(":: Reference and query before construction");
            LOG.debug(String.format(":: read=%s", sourceQuery));
        }

        applyCigar();
        applyMd();
        clipRefAndQuery();
        alignedLength = query.length();
        queryAlignedLength = alignedLength - numDeletions;
        targetAlignedLength = alignedLength - numInsertions;
        score = alignedLength - numDeletions - numInsertions - numMisMatches;

        // Figure out start of alignment and alignment length, by observing mismatches at head and tail
        if (query.length() != ref.length()) {
            LOG.error("ERROR! reconstructed reads and refs don't match in size!!");
        }
    }

    protected void applyCigar() {
        if (debug && LOG.isDebugEnabled()) {
            LOG.debug(String.format(":: Applying cigar=%s", cigar));
        }
        int posInReads = 0;
        numInsertions = 0;
        numDeletions = 0;
        numMisMatches = 0;
        Matcher matcher = CIGAR_REGEX.matcher(cigar);
        while (matcher.find()) {
            final int length = Integer.parseInt(matcher.group(1));
            final char op = matcher.group(2).charAt(0);
            switch (op) {
                case 'S':
                    // Soft clipping
                    for (int i = 0; i < length; i++) {
                        ref.append('-');
                        query.append('-');
                        qual.append((char) minQualValue); // min quality
                    }
                    //  posInReads += length;
                    break;
                case 'M':
                    // Account for matches AND mismatches. Any mis-matches will be fixed in applyMd()

                    query.append(sourceQuery.substring(posInReads, posInReads + length));
                    if (sourceQual.length() != 0) {
                        qual.append(sourceQual.substring(posInReads, posInReads + length));
                    }
                    ref.append(sourceQuery.substring(posInReads, posInReads + length));
                    posInReads += length;
                    break;
                case 'I':
                    query.append(sourceQuery.substring(posInReads, posInReads + length));
                    if (sourceQual.length() != 0) {
                        qual.append(sourceQual.substring(posInReads, posInReads + length));
                    }
                    for (int i = 0; i < length; i++) {
                        ref.append('-');
                    }
                    numInsertions += length;
                    posInReads += length;
                    break;
                case 'D':
                    for (int i = 0; i < length; i++) {
                        query.append('-');
                        if (sourceQual.length() != 0) {
                            // Minimum qual, placing min value here but it shouldn't get written to
                            // sequence variations
                            qual.append((char) minQualValue);
                        }
                        ref.append('?');
                    }
                    numDeletions += length;
                    break;
            }
        }
        for (int i = 0; i < query.length(); i++) {
            if (query.charAt(i) == '-' && ref.charAt(i) == '-') {
                numLeftClipped++;
            } else {
                break;
            }
        }
        for (int i = query.length() - 1; i >= 0; i--) {
            if (query.charAt(i) == '-' && ref.charAt(i) == '-') {
                numRightClipped++;
            } else {
                break;
            }
        }
    }

    private void applyMd() {
        // My RE is simplified from the SAM spec but performs the same task
        // the major difference being mine would allow 5ACG where the current
        // spec would require 5A0C0G0 (which mine will still work with just fine).
        if (debug && LOG.isDebugEnabled()) {
            LOG.debug(String.format(":: Applying md=%s", md));
        }
        int position =  numLeftClipped;
        final Matcher matcher = MD_REGEX.matcher(md);
        while (matcher.find()) {
            final String mdPart = matcher.group();
            if (NUMERIC_REGEX.matcher(mdPart).matches()) {

                final int length = Integer.parseInt(mdPart);
                position += length;

            } else if (mdPart.charAt(0) == '^') {
                if (ref.charAt(position) == '-') {
                    position++;
                }
                // Adjust the ref with these characters, ignoring the ^ character so start at 1
                for (int i = 1; i < mdPart.length(); i++) {
                    ref.setCharAt(position++, mdPart.charAt(i));
                }
            } else {
                // The regex should only allow a single character here, but we'll accept multiple
                for (int i = 0; i < mdPart.length(); i++) {
                    if (ref.charAt(position) == '-') {
                        position++;
                    }
                    ref.setCharAt(position++, Character.toLowerCase(mdPart.charAt(i)));
                    numMisMatches++;
                }
            }
        }
    }

    private void clipRefAndQuery() {
        if (numLeftClipped > 0 || numRightClipped > 0) {
            LOG.debug(":: Reference and query pre-clipping");
            debugSequences(false);
            if (numRightClipped > 0) {
                query.setLength(query.length() - numRightClipped);
                ref.setLength(ref.length() - numRightClipped);
                if (qual.length() > 0) {
                    qual.setLength(qual.length() - numRightClipped);
                }
            }
            if (numLeftClipped > 0) {
                queryPosition += numLeftClipped;
                query.delete(0, numLeftClipped);
                ref.delete(0, numLeftClipped);
                if (qual.length() > 0) {
                    qual.delete(0, numLeftClipped);
                }
            }
        }
    }

    private void findSequenceVariations() {
        char refChar, queryChar;
        int readIndex, refPosition;
        int genomicLength = ref.length();
        int paddedLength = numLeftClipped + genomicLength + numRightClipped;
        boolean tooBig = false;
        int tooBigReadIndex = 0, tooBigRefPosition = 0;
        refPositions.size(10);

        readIndexes.size(paddedLength);
        refPositions.size(paddedLength);

        refPosition = 0;
        readIndex = 0;
        for (int i = 0; i < paddedLength; i++) {
            if (i < numLeftClipped) {
                // In alignment padding refPosition doesn't increment
                // but read position does
                readIndex++;
            } else if (i >= (genomicLength + numLeftClipped)) {
                // In alignment padding/clipping refPosition doesn't increment
                // but read position does
                readIndex++;
            } else {
                refChar = Character.toUpperCase(ref.charAt(i - numLeftClipped));
                if (refChar != '-') {
                    refPosition++;
                }
                if (reverseStrand) {
                    final int index = genomicLength - (i - numLeftClipped) - 1;
                    queryChar = Character.toUpperCase(query.charAt(index));
                } else {
                    final int index = i - numLeftClipped;
                    queryChar = Character.toUpperCase(query.charAt(index));
                }
                if (queryChar != '-') {
                    readIndex++;
                }
            }
            refPositions.set(i, refPosition);
            if (reverseStrand) {
                readIndexes.set(paddedLength - i - 1, readIndex);
            } else {
                readIndexes.set(i, readIndex);
            }
        }

        if (debug && LOG.isDebugEnabled()) {
            debugSequences();
            logval.setLength(0);
            logval.append("::  pos=");
            for (int i = 0; i < paddedLength; i++) {
                logval.append(String.format("%d", refPositions.get(i) % 10));
            }
            LOG.debug(logval.toString());
            logval.setLength(0);
            logval.append("::   ri=");
            for (int i = 0; i < paddedLength; i++) {
                logval.append(String.format("%d", readIndexes.get(i) % 10));
            }
            LOG.debug(logval.toString());

            LOG.debug("ref with positions");
            logval.setLength(0);
            for (int i = 0; i < paddedLength; i++) {
                if (i < numLeftClipped) {
                    refChar = '_';
                } else if (i >= (genomicLength + numLeftClipped)) {
                    refChar = '_';
                } else {
                    refChar = ref.charAt(i - numLeftClipped);
                }
                logval.append(String.format("%02d:%c:%02d  ", i, refChar, refPositions.get(i)));
            }
            LOG.debug(logval.toString());

            LOG.debug("read with positions");
            logval.setLength(0);
            for (int i = 0; i < paddedLength; i++) {
                if (i < numLeftClipped) {
                    queryChar = '_';
                } else if (i >= (numLeftClipped + genomicLength)) {
                    queryChar = '_';
                } else {
                    queryChar = query.charAt(i - numLeftClipped);
                }
                logval.append(String.format("%02d:%c:%02d  ", i, queryChar, readIndexes.get(i)));
            }
            LOG.debug(logval.toString());
            LOG.debug(String.format("numLeftClipped=%d, numRightClipped=%d", numLeftClipped, numRightClipped));
        }

        for (int queryI = numLeftClipped; queryI < numLeftClipped + genomicLength; queryI++) {
            refPosition = refPositions.get(queryI);
            readIndex = readIndexes.get(queryI);
            final int i = queryI - numLeftClipped;
            if (readIndex > queryLength && !tooBig) {
                tooBig = true;
                tooBigReadIndex = readIndex;
                tooBigRefPosition = refPosition;
            }
            refChar = Character.toUpperCase(ref.charAt(i));
            queryChar = Character.toUpperCase(query.charAt(i));
            boolean hasQual;
            byte qualChar;
            // We check queryChar != '-' because we don't have a qual score on deletions
            if (qual.length() > 0 && queryChar != '-') {
                hasQual = true;
                qualChar = qualityEncoding.asciiEncodingToPhredQualityScore(qual.charAt(i));
            } else {
                hasQual = false;
                qualChar = minQualValue;
            }
            if (refChar != queryChar) {
                sequenceVariations.add(new SamSequenceVariation(refPosition, refChar, readIndex, queryChar, hasQual, qualChar));
            }
        }

        if (tooBig) {
            if (debug && LOG.isDebugEnabled()) {
                LOG.debug(String.format(" *** readIndex [%d] or refPosition [%d] is too large! ***",
                        tooBigReadIndex, tooBigRefPosition));
                LOG.debug(String.format(">%d", queryIndex));
                if (sourceQuery.length() > 0) {
                    LOG.debug(String.format("%s", sourceQuery));
                } else {
                    LOG.debug("sourceQuery was NULL");
                }
            }
        }
    }

    private void debugSequences() {
        debugSequences(true);
    }

    private void debugSequences(final boolean printClipChars) {
        if (!debug) {
            return;
        }
        LOG.debug(String.format(":: paddingLeft=%d, paddingRight=%d", numLeftClipped, numRightClipped));
        logval.setLength(0);
        if (qual.length() > 0) {
            logval.append(":: qual=");
            if (printClipChars) {
                for (int i = 0; i < numLeftClipped; i++) {
                    logval.append("_");
                }
            }
            for (int i = 0; i < qual.length(); i++) {
                if (i > 0) {
                    logval.append(":");
                }
                logval.append(String.format("%d", (int) qual.charAt(i)));
            }
            if (printClipChars) {
                for (int i = 0; i < numRightClipped; i++) {
                    logval.append("_");
                }
            }
        } else {
            logval.append(":: qual=none");
        }
        LOG.debug(logval);

        logval.setLength(0);
        logval.append(":: ref =");
        if (printClipChars) {
            for (int i = 0; i < numLeftClipped; i++) {
                logval.append("_");
            }
        }
        logval.append(ref.toString());
        if (printClipChars) {
            for (int i = 0; i < numRightClipped; i++) {
                logval.append("_");
            }
        }
        LOG.debug(logval.toString());

        logval.setLength(0);
        logval.append(":: read=");
        if (printClipChars) {
            for (int i = 0; i < numLeftClipped; i++) {
                logval.append("_");
            }
        }
        logval.append(query.toString());
        if (printClipChars) {
            for (int i = 0; i < numRightClipped; i++) {
                logval.append("_");
            }
        }
        LOG.debug(logval.toString());
    }

    /**
     * Get the quality encoding scale used for the input fastq file.
     *
     * @return the quality encoding scale used for the input fastq file
     */
    public QualityEncoding getQualityEncoding() {
        return qualityEncoding;
    }

    /**
     * Set the quality encoding scale to be used for the input fastq file.
     * Acceptable values are "Illumina", "Sanger", and "Solexa".
     *
     * @param qualityEncoding the quality encoding scale to be used for the input fastq file
     */
    public void setQualityEncoding(final QualityEncoding qualityEncoding) {
        this.qualityEncoding = qualityEncoding;
    }

    public void setQueryPosition(final int queryPosition) {
        this.queryPosition = queryPosition;
    }

    /**
     * When appending two cigar strings, if the left ends with the same code as the right stars with,
     * this will combine those two into a single code. Such as appendTo="2S45M" + appendFrom="25M2S"
     * the result will be appendTo=="2S70M2S".
     *
     * @param appendTo   the cigar string we are appending to
     * @param appendFrom the cigar string we are appending from
     */
    public static void appendCigar(final MutableString appendTo, final MutableString appendFrom) {
        // A couple simple cases
        if (appendFrom.length() == 0 || appendTo.length() == 0) {
            appendTo.append(appendFrom);
            return;
        }

        String appendToMatchNum = null;
        char appendToMatchCode = '\0';
        final Matcher appendToMatcher = LAST_CIGAR_PATTERN.matcher(appendTo);
        if (appendToMatcher.find()) {
            appendToMatchNum = appendToMatcher.group(1);
            appendToMatchCode = appendToMatcher.group(2).charAt(0);
        }
        String appendFromMatchNum = null;
        char appendFromMatchCode = '\0';
        final Matcher appendFromMatcher = FIRST_CIGAR_PATTERN.matcher(appendFrom);
        if (appendFromMatcher.find()) {
            appendFromMatchNum = appendFromMatcher.group(1);
            appendFromMatchCode = appendFromMatcher.group(2).charAt(0);
        }
        if (appendToMatchNum == null || appendFromMatchNum == null || appendToMatchCode != appendFromMatchCode) {
            // Either appendTo doesn't end in a number of appendFrom doesn't start with a number
            appendTo.append(appendFrom);
        } else {
            final int appendToVal = Integer.parseInt(appendToMatchNum);
            final int appendFromVal = Integer.parseInt(appendFromMatchNum);
            appendTo.setLength(appendTo.length() - appendToMatchNum.length() - 1);
            appendTo.append(appendToVal + appendFromVal).append(appendFromMatchCode);
            appendTo.append(appendFrom.substring(appendFromMatchNum.length() + 1));
        }
    }

    /**
     * When one has two MD:Z mismatch strings, if the left ends with a number and the right ends with a number,
     * you cannot just concatenate left+right, you must take the left number and add it to the right number.
     * If appendTo="10TCC25" + appendFrom="15^TCC10" will result in appendTo=="10TCC40^TCC10".
     *
     * @param appendTo   the mismatchString we are appending to
     * @param appendFrom the mismatchString we are appending from
     */
    public static void appendMismatches(final MutableString appendTo, final MutableString appendFrom) {
        // A couple simple cases
        if (appendFrom.length() == 0 || appendTo.length() == 0) {
            appendTo.append(appendFrom);
            return;
        }

        String appendToMatch = null;
        final Matcher appendToMatcher = LAST_NUMBER_PATTERN.matcher(appendTo);
        if (appendToMatcher.find()) {
            appendToMatch = appendToMatcher.group(0);
        }
        String appendFromMatch = null;
        final Matcher appendFromMatcher = FIRST_NUMBER_PATTERN.matcher(appendFrom);
        if (appendFromMatcher.find()) {
            appendFromMatch = appendFromMatcher.group(0);
        }
        if (appendToMatch == null || appendFromMatch == null) {
            // Either appendTo doesn't end in a number of appendFrom doesn't start with a number
            appendTo.append(appendFrom);
        } else {
            final int appendToVal = Integer.parseInt(appendToMatch);
            final int appendFromVal = Integer.parseInt(appendFromMatch);
            appendTo.setLength(appendTo.length() - appendToMatch.length());
            appendTo.append(appendToVal + appendFromVal);
            appendTo.append(appendFrom.substring(appendFromMatch.length()));
        }
    }
}
