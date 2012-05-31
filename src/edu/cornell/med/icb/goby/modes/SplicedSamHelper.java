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

package edu.cornell.med.icb.goby.modes;

import edu.cornell.med.icb.goby.reads.QualityEncoding;
import edu.cornell.med.icb.goby.util.pool.NullResettableObjectPool;
import edu.cornell.med.icb.goby.util.pool.ResettableObjectPoolInterface;
import it.unimi.dsi.Util;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A facade over SamHelper that handles spliced CIGAR strings and generates multiple alignment entries. This is a draft
 * that needs to be better tested, but illustrates how to extend SamHelper functionality to support spliced BAM input.
 *
 * @author Fabien Campagne
 *         Date: 3/29/12
 *         Time: 10:51 AM
 */
public class SplicedSamHelper {
    private static final Logger LOG = Logger.getLogger(SplicedSamHelper.class);

    private final boolean debug;

    private int numEntries;
    private int cursorIndex;
    private ObjectArrayList<SamHelper> helpers = new ObjectArrayList<SamHelper>();
    private final ResettableObjectPoolInterface<SamHelper> helperPool = new NullResettableObjectPool<SamHelper>() {
        @Override
        public SamHelper makeObject() {
            return new SamHelper();
            // The helper factory was doing this, but the helper factory was no longer being used.
            // if (encoding != null) {
            //     samHelper.setQualityEncoding(encoding);
            // }
        }
    };
    private static final Pattern CIGAR_REGEX = Pattern.compile("([0-9]+)([SMIDN])");
    private static final Pattern MD_REGEX = Pattern.compile("([0-9]+|[ACGTN]|\\^[ACGTN]+)");
    private static final Pattern NUMERIC_REGEX = Pattern.compile("^[0-9]+$");
    private boolean usingGenome;
    private int queryLength;

    private QualityEncoding encoding;


    /**
     * Compute limits from cigar and mdString. Each N stretch in the cigar string separates two alignment segments.
     * The Limits instances describe individual segments.
     *
     * @param position the position
     * @param cigar the cigar string
     * @param mdString the md string
     * @return the limits
     */
    protected Limits[] getLimits(final int position, final String cigar, final CharSequence mdString) {
        if (debug && LOG.isDebugEnabled()) {
            LOG.debug(String.format(":: Applying cigar=%s", cigar));
        }
        final ObjectArrayList<Limits> list = new ObjectArrayList<Limits>();

        analyzeCigar(position, cigar, list);
        if (mdString != null) {
            analyzeMd(position, mdString.toString(), list);
        }
        return list.toArray(new Limits[list.size()]);
    }

    private void analyzeMd(final int position, final String mdString, final ObjectArrayList<Limits> list) {
        for (final Limits limit : list) {
            adjust(mdString, limit);
        }
        if (list.size() == 1) {
            list.get(0).md = mdString;
        }
    }

    private void adjust(final String mdString, final Limits limit) {
        final Matcher matcher = MD_REGEX.matcher(mdString);
        int mdIndex = 0;
        int previousMdIndex = 0;
        int positionInRead = 0;
        int previousPositionInRead = 0;
        while (matcher.find()) {
            final String mdPart = matcher.group();
            if (NUMERIC_REGEX.matcher(mdPart).matches()) {
                final int baseLength = Integer.parseInt(mdPart);
                positionInRead += baseLength;
                mdIndex += mdPart.length();
                if (overlaps(limit, previousPositionInRead, positionInRead)) {
                    final int matchingSpanInLimit = Math.min(limit.readEnd - limit.readStart, baseLength);
                    limit.md = limit.md + Integer.toString(matchingSpanInLimit);
                }
                previousPositionInRead = positionInRead;
            } else {
                mdIndex += mdPart.length();
                positionInRead += mdPart.length();
                if (overlaps(limit, previousPositionInRead, positionInRead)) {

                    limit.md = mdString.substring(previousMdIndex, mdIndex);

                }
                previousPositionInRead = positionInRead;
            }
            previousMdIndex = mdIndex;
        }

    }

    private boolean overlaps(final Limits limit, final int limitStartA, final int limitEndA) {
        final int limitStartB = limit.readStart;
        final int limitEndB = limit.readEnd;
        return limitStartA >= limitStartB && limitEndA <= limitEndB;
    }

    private void analyzeCigar(int position, final String cigar, final ObjectArrayList<Limits> list) {
        final Matcher matcher = CIGAR_REGEX.matcher(cigar);
        int previousCigarIndex = 0;
        int cigarIndex = 0;
        int previousPosition = position;
        int previousPositionInRead = 0;
        int positionInRead = 0;
        int trim = 0;
        while (matcher.find()) {
            final int cigarLength = matcher.group(1).length() + matcher.group(2).length();
            final int readBasesLength = Integer.parseInt(matcher.group(1));
            final char op = matcher.group(2).charAt(0);

            switch (op) {

                case 'N':
                    //  System.out.println("new CIGAR: "+cigar.substring(previousCigarIndex, cigarIndex));
                    final Limits limits = new Limits(previousPosition, previousCigarIndex, cigarIndex, previousPositionInRead, positionInRead, previousPosition, position);
                    limits.setTrim(trim);
                    trim = 0;
                    list.add(limits);
                    previousCigarIndex = cigarIndex + cigarLength; // we exclude the N group from the limits.
                    previousPositionInRead = positionInRead;
                    position += readBasesLength;
                    previousPosition = position;
                    break;

                case 'S':

                    positionInRead += readBasesLength;
                    trim = readBasesLength;
                    //   previousPositionInRead += readBasesLength;
                    break;
                default:

                    positionInRead += readBasesLength;
                    position += readBasesLength;
                    break;
            }
            cigarIndex += cigarLength;

        }
        final Limits newLimit = new Limits(previousPosition, previousCigarIndex, cigarIndex, previousPositionInRead,
                positionInRead, previousPosition, position);
        newLimit.setTrim(trim);
        list.add(newLimit);
    }

    public void setQualityEncoding(final QualityEncoding qualityEncoding) {
        encoding = qualityEncoding;
    }

    public SplicedSamHelper() {
        // don't even dare go through the debugging code if log4j was not configured. The debug code
        // is way too slow to run unintentionally in production!
        debug = Util.log4JIsConfigured();
        helpers = new ObjectArrayList<SamHelper>();
        reset();

    }


    public final void reset() {
        try {
            for (final SamHelper h : helpers) {
                helperPool.returnObject(h);
            }
        } catch (Exception e) {
            LOG.error("Unable to return object to pool ", e);
        }
        cursorIndex = 0;
        numEntries = 1;
        helpers.clear();
        refSequence.setLength(0);
        usingGenome = false;
    }

    private static final String N_STRING = "N";

    public void setSourceWithReference(final int queryIndex, final SAMRecord samRecord, final String sourceReference) {
        final String cigarString = samRecord.getCigarString();
        final int position = samRecord.getAlignmentStart();    // one-based
        refSequence.append(sourceReference);
        usingGenome = true;
        final Limits[] limits = getLimits(position, cigarString, null);
        final CharSequence sourceQuery = samRecord.getReadString();
        final CharSequence sourceQual = samRecord.getBaseQualityString();
        final boolean reverseStrand = samRecord.getReadNegativeStrandFlag();
        numEntries = limits.length;
        initializeHelpers();
        queryLength = samRecord.getReadLength();
        for (int i = 0; i < numEntries; i++) {
            final Limits limit = limits[i];
            final int refStartIndex = limit.refStart - position;
            final int refEndIndex = refStartIndex + limit.refEnd - limit.refStart;
            try {
                final CharSequence sourceRef = refSequence.subSequence(refStartIndex, Math.min(refEndIndex, sourceReference.length() - 1));
                helpers.get(i).setSourceWithReference(queryIndex,
                        sourceRef,
                        sourceQuery.subSequence(limit.readStart, limit.readEnd),
                        sourceQual.subSequence(limit.readStart, limit.readEnd),
                        limit.position,
                        reverseStrand
                );
                final int queryPosition = (i != numEntries - 1 ? limit.trim : 0) + limit.readStart;
                helpers.get(i).setQueryPosition(queryPosition);

            } catch (IndexOutOfBoundsException e) {
                System.out.printf("Another exception: refStartIndex=%d refEndIndex=%d refLength=%d  cigar=%s %s %n",
                        refStartIndex, refEndIndex, sourceReference.length(), samRecord.getCigarString(), e);
            }
        }
    }

    MutableString refSequence = new MutableString();

    /**
     * Set the source
     * @param queryIndex query index
     * @param sourceQuery source query
     * @param sourceQual source qualities
     * @param cigar cigar string
     * @param md md:z string
     * @param position      one-based position
     * @param reverseStrand if reverse strand
     * @param readLength the read length
     */
    public void setSource(final int queryIndex, final CharSequence sourceQuery, final CharSequence sourceQual,
                          final CharSequence cigar, final CharSequence md, final int position,
                          final boolean reverseStrand,
                          final int readLength) {

        usingGenome = false;
        final String cigarString = cigar.toString();
        queryLength = readLength;

        if (!cigarString.contains(N_STRING)) {
            numEntries = 1;
            cursorIndex = 0;
            initializeHelpers();
            helpers.get(0).setSource(queryIndex, sourceQuery, sourceQual, cigar, md, position, reverseStrand, queryLength);
        } else {

            final Limits[] limits = getLimits(position, cigarString, md);
            numEntries = limits.length;
            initializeHelpers();
            for (int i = 0; i < numEntries; i++) {
                final Limits limit = limits[i];
                helpers.get(i).setSource(queryIndex,
                        sourceQuery.subSequence(limit.readStart, limit.readEnd),
                        sourceQual.subSequence(limit.readStart, limit.readEnd),
                        cigarString.substring(limit.cigarStart, limit.cigarEnd),
                        limit.md,
                        limit.position,
                        reverseStrand, queryLength
                );
                final int queryPosition = (i != numEntries - 1 ? limit.trim : 0) + limit.readStart;
                helpers.get(i).setQueryPosition(queryPosition);
            }
        }
    }

    private void initializeHelpers() {
        try {
            // setup the helpers instances, getting them from the pool as needed:
            helpers.clear();
            for (int j = 0; j < numEntries; j++) {

                helpers.add(helperPool.borrowObject());
            }
        } catch (Exception e) {
            LOG.error("Unable to borrow object from pool ", e);
        }
    }

    protected static class Limits {
        public int readStart;
        public int readEnd;
        public int cigarStart;
        public int cigarEnd;
        public int mdStart = Integer.MAX_VALUE;
        public int mdEnd = Integer.MIN_VALUE;
        public String md = "";
        public int refStart;
        public int refEnd;
        private int trim;

        private Limits(final int position, final int cigarStart, final int cigarEnd, final int readStart,
                       final int readEnd, final int refStart, final int refEnd) {
            this.cigarStart = cigarStart;
            this.cigarEnd = cigarEnd;
            this.readStart = readStart;
            this.readEnd = readEnd;
            this.refStart = refStart;
            this.refEnd = refEnd;
            this.position = position;
        }

        /**
         * one-based position.
         */
        public int position;

        public void setTrim(final int trim) {
            this.trim = trim;
        }
    }


    /**
     * Return the number of Goby alignment entries represented by the source. When the NM field contains Ns, a SAM entry
     * may be split into two or more Goby alignment entries linked via SpliceLinks. This method returns the number of
     * such entries. See the setEntryCursor(index) method to choose which result entry to get data about.
     *
     * @return the maximum cursor index+1.
     */
    public int getNumEntries() {
        return numEntries;
    }

    /**
     * Choose which alignment entry to return data about.
     *
     * @param cursorIndex an index less than getNumEntries() and at least zero.
     */
    public void setEntryCursor(final int cursorIndex) {
        this.cursorIndex = cursorIndex;
    }

    public int getAlignedLength() {
        return helpers.get(cursorIndex).getAlignedLength();
    }

    public MutableString getCigar() {
        return helpers.get(cursorIndex).getCigar();
    }

    public MutableString getMd() {
        return helpers.get(cursorIndex).getMd();
    }

    public byte getMinQualValue() {
        return helpers.get(cursorIndex).getMinQualValue();
    }

    public int getNumDeletions() {
        return helpers.get(cursorIndex).getNumDeletions();
    }

    public int getNumInsertions() {
        return helpers.get(cursorIndex).getNumInsertions();
    }

    public int getNumLeftClipped() {
        return helpers.get(cursorIndex).getNumLeftClipped();
    }

    public int getNumMisMatches() {
        return helpers.get(cursorIndex).getNumMisMatches();
    }

    public int getNumRightClipped() {
        return helpers.get(cursorIndex).getNumRightClipped();
    }


    public byte[] getSourceQualAsBytes() {
        if (cursorIndex==0) {
        return helpers.get(cursorIndex).getSourceQualAsBytes();
        } else {
            return null;
        }
    }

    /**
     * Return zero-based position.
     *
     * @return zero-based position.
     */
    public int getPosition() {
        return helpers.get(cursorIndex).getPosition();
    }

    public MutableString getQual() {

        return helpers.get(cursorIndex).getQual();

    }

    public QualityEncoding getQualityEncoding() {
        return helpers.get(cursorIndex).getQualityEncoding();
    }

    public MutableString getQuery() {
        return helpers.get(cursorIndex).getQuery();
    }

    public int getQueryAlignedLength() {
        return helpers.get(cursorIndex).getQueryAlignedLength();
    }

    public int getQueryIndex() {
        return helpers.get(cursorIndex).getQueryIndex();
    }

    public int getQueryLength() {
        return queryLength;
    }

    public int getQueryPosition() {
        return helpers.get(cursorIndex).getQueryPosition();
    }

    public MutableString getRef() {
        return helpers.get(cursorIndex).getRef();
    }

    public int getScore() {
        return helpers.get(cursorIndex).getScore();
    }

    public List<SamSequenceVariation> getSequenceVariations() {
        return helpers.get(cursorIndex).getSequenceVariations();
    }

    public MutableString getSourceQual() {
        if (cursorIndex==0) {
        return helpers.get(cursorIndex).getSourceQual();
        } else {
            return null;
        }
    }

    public MutableString getSourceQuery() {
        return helpers.get(cursorIndex).getSourceQuery();
    }

    public int getTargetAlignedLength() {
        return helpers.get(cursorIndex).getTargetAlignedLength();
    }

    public boolean isReverseStrand() {
        return helpers.get(cursorIndex).isReverseStrand();
    }


}
