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

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.modes.SamHelper;
import it.unimi.dsi.lang.MutableString;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.apache.commons.lang.ArrayUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Class for comparing two SAMRecords. Used with both tests and
 * with SAMComparisonMode.
 */
public class SamComparison {

    /**
     * Configuration values.
     */
    private boolean mappedQualitiesPreserved;
    private boolean softClipsPreserved;
    private boolean checkMate;
    private boolean canonicalMdzForComparison;

    /**
     * Running totals, etc.
     */
    private int readNum;
    private int comparisonFailureCount;

    /**
     * Data to compare.
     */
    private SAMRecord expectedSamRecord;
    private SAMRecord actualSamRecord;
    private Alignments.AlignmentEntry gobyAlignment;

    /**
     * State.
     */
    private final MutableString comparisonErrorDump;
    private List<String> comparisonFailures;
    private boolean initialized;

    /**
     * Create a SamComparison.
     */
    public SamComparison() {
        comparisonErrorDump = new MutableString();
        comparisonFailures = new ArrayList<String>();
        initialized = false;
        mappedQualitiesPreserved = false;
        softClipsPreserved = false;
        checkMate = false;
        canonicalMdzForComparison = true;
    }

    /**
     * Call this one time before at the start of a large number of comparisons (such as before comparing
     * one entire file with another). Calling the first time is optional as it will be called if a comparison
     * is done and it has never been called.
     */
    public void reset() {
        readNum = 0;
        comparisonFailureCount = 0;
    }

    /**
     * Get if it is assumed that the compact file created from the BAM/SAM
     * file preserved mapped qualities.
     * @return if it is assumed ...
     */
    public boolean isMappedQualitiesPreserved() {
        return mappedQualitiesPreserved;
    }

    /**
     * Set if it is assumed that the compact file created from the BAM/SAM
     * file preserved mapped qualities.
     * @param mappedQualitiesPreserved if it is assumed...
     */
    public void setMappedQualitiesPreserved(final boolean mappedQualitiesPreserved) {
        this.mappedQualitiesPreserved = mappedQualitiesPreserved;
    }

    /**
     * Get if it is assumed that the compact file created from the BAM/SAM
     * file preserved soft clips.
     * @return if it is assumed ...
     */
    public boolean isSoftClipsPreserved() {
        return softClipsPreserved;
    }

    /**
     * Set if it is assumed that the compact file created from the BAM/SAM
     * file preserved soft clips.
     * @param softClipsPreserved if it is assumed ...
     */
    public void setSoftClipsPreserved(final boolean softClipsPreserved) {
        this.softClipsPreserved = softClipsPreserved;
    }

    /**
     * Get if the details about mate reads will be checked.
     * If the source SAM/BAM file is a complete file you can set this to true,
     * if you are using an incomplete source SAM/BAM file, this should be
     * set to false. Default is false.
     * @return if mates will be checked
     */
    public boolean isCheckMate() {
        return checkMate;
    }

    /**
     * Set if the details about mate reads will be checked.
     * If the source SAM/BAM file is a complete file you can set this to true,
     * if you are using an incomplete source SAM/BAM file, this should be
     * set to false. Default is false.
     * @return if mates will be checked
     */
    public void setCheckMate(final boolean checkMate) {
        this.checkMate = checkMate;
    }

    /**
     * Get if canonical MD:Z comparisons will be made.
     * When true, the source and destination MD:Z values will be passed through an algorithm
     * to make them canonical (place 0's in places where 0's should exist but might not).
     * By default this is enabled.
     * @return if ...
     */
    public boolean isCanonicalMdzForComparison() {
        return canonicalMdzForComparison;
    }

    /**
     * Set if canonical MD:Z comparisons will be made.
     * When true, the source and destination MD:Z values will be passed through an algorithm
     * to make them canonical (place 0's in places where 0's should exist but might not).
     * By default this is enabled.
     * @param canonicalMdzForComparison if ...
     */
    public void setCanonicalMdzForComparison(final boolean canonicalMdzForComparison) {
        this.canonicalMdzForComparison = canonicalMdzForComparison;
    }

    /**
     * Return how many reads have been compared since reset() was last called.
     * @return how many...
     */
    public int getReadNum() {
        return readNum;
    }

    /**
     * Return how many comparison failures have been found since reset() was last called.
     * @return how many...
     */
    public int getComparisonFailureCount() {
        return comparisonFailureCount;
    }

    /**
     * Get the expected sam record. This contains the values we EXPECT TO FIND, ie, the source SAM/BAM file.
     * @return the expected sam record.
     */
    public SAMRecord getExpectedSamRecord() {
        return expectedSamRecord;
    }

    /**
     * Set the expected sam record. This contains the values we EXPECT TO FIND, ie, the source SAM/BAM file.
     * @param expectedSamRecord the expected sam record.
     */
    public void setExpectedSamRecord(final SAMRecord expectedSamRecord) {
        this.expectedSamRecord = expectedSamRecord;
    }

    /**
     * Get the ACTUAL sam record. This contains the values we will compare against the expected, ie, the destination
     * SAM/BAM file.
     * @return the actual sam record.
     */
    public SAMRecord getActualSamRecord() {
        return actualSamRecord;
    }

    /**
     * Set the ACTUAL sam record. This contains the values we will compare against the expected, ie, the destination
     * SAM/BAM file.
     * @param actualSamRecord the actual sam record.
     */
    public void setActualSamRecord(final SAMRecord actualSamRecord) {
        this.actualSamRecord = actualSamRecord;
    }

    /**
     * Get the Goby Compact Alignment (optional). This should be the intermediate record that was created
     * from the expected record and was used to create the actual record. Providing this will generate some additional
     * output when there are comparison failures.
     * @return the goby alignment
     */
    public Alignments.AlignmentEntry getGobyAlignment() {
        return gobyAlignment;
    }

    /**
     * Set the Goby Compact Alignment (optional). This should be the intermediate record that was created
     * from the expected record and was used to create the actual record. Providing this will generate some additional
     * output when there are comparison failures.
     * @param gobyAlignment the goby alignment
     */
    public void setGobyAlignment(final Alignments.AlignmentEntry gobyAlignment) {
        this.gobyAlignment = gobyAlignment;
    }

    /**
     * Compare expectedSamRecord vs actualSamRecord. Output details if differences are found.
     * @return true if the two records are found to be the same
     */
    public boolean compareSamRecords() {
        if (!initialized) {
            initialized = true;
            reset();
        }
        comparisonFailures.clear();
        compareField("Positions don't match", expectedSamRecord.getAlignmentStart(), actualSamRecord.getAlignmentStart());
        compareField("Ref Index don't match", expectedSamRecord.getReferenceIndex(), actualSamRecord.getReferenceIndex());
        compareField("Flags don't match", expectedSamRecord.getFlags(), actualSamRecord.getFlags());
        compareField("Mapping quality doesn't match", expectedSamRecord.getMappingQuality(), actualSamRecord.getMappingQuality());
        compareField("Read paired flag doesn't match", expectedSamRecord.getReadPairedFlag(), actualSamRecord.getReadPairedFlag());
        compareField("Read length doesn't match", expectedSamRecord.getReadLength(), actualSamRecord.getReadLength());
        if (checkMate) {
            if (expectedSamRecord.getReadPairedFlag()) {
                compareField("Read mate unmapped doesn't match", expectedSamRecord.getMateUnmappedFlag(), actualSamRecord.getMateUnmappedFlag());
                if (!expectedSamRecord.getMateUnmappedFlag()) {
                    compareField("Mate alignment start doesn't match", expectedSamRecord.getMateAlignmentStart(), actualSamRecord.getMateAlignmentStart());
                    compareField("Mate alignment index doesn't match", expectedSamRecord.getMateReferenceIndex(), actualSamRecord.getMateReferenceIndex());
                    compareField("Inferred insert size doesn't match", expectedSamRecord.getInferredInsertSize(), actualSamRecord.getInferredInsertSize());
                }
            }
        }
        compareField("Positive/negative strand doesn't match", expectedSamRecord.getReadNegativeStrandFlag(), actualSamRecord.getReadNegativeStrandFlag());
        compareField("Cigars don't match", expectedSamRecord.getCigarString(), actualSamRecord.getCigarString());

        final String eMdz;
        final String aMdz;
        if (canonicalMdzForComparison) {
            eMdz = SamHelper.canonicalMdz(expectedSamRecord.getStringAttribute("MD"));
            aMdz = SamHelper.canonicalMdz(actualSamRecord.getStringAttribute("MD"));
        } else {
            eMdz = expectedSamRecord.getStringAttribute("MD");
            aMdz = actualSamRecord.getStringAttribute("MD");
        }
        compareField("MD:Z doesn't match", eMdz, aMdz);

        final String eRead;
        final String aRead;
        if (softClipsPreserved) {
            eRead = expectedSamRecord.getReadString();
            aRead = actualSamRecord.getReadString();
        } else {
            int clipLeft = 0;
            int clipRight = 0;
            final List<CigarElement> eCigarElements = expectedSamRecord.getCigar().getCigarElements();
            if (!eCigarElements.isEmpty()) {
                final CigarElement firstCigar = eCigarElements.get(0);
                if (firstCigar.getOperator() == CigarOperator.S) {
                    clipLeft = firstCigar.getLength();
                }
            }
            if (eCigarElements.size() > 1) {
                final CigarElement lastCigar = eCigarElements.get(eCigarElements.size() - 1);
                if (lastCigar.getOperator() == CigarOperator.S) {
                    clipRight = lastCigar.getLength();
                }
            }
            eRead = expectedSamRecord.getReadString().substring(clipLeft, expectedSamRecord.getReadLength() - clipRight);
            aRead = actualSamRecord.getReadString().substring(clipLeft, actualSamRecord.getReadLength() - clipRight);
            compareField("Clipped size doesn't compute",
                    actualSamRecord.getReadLength(), aRead.length() + clipLeft + clipRight);
        }

        compareField("Reads didn't match", eRead, aRead);
        if (mappedQualitiesPreserved) {
            compareField("Quality didn't match", expectedSamRecord.getBaseQualityString(), actualSamRecord.getBaseQualityString());
        } else {
            if (gobyAlignment != null) {
                final int readLength = expectedSamRecord.getReadLength();
                for (final Alignments.SequenceVariation seqvar : gobyAlignment.getSequenceVariationsList()) {
                    final String to = seqvar.getTo();
                    int i = 0;
                    for (final char toChar : to.toCharArray()) {
                        if (toChar != '-') {
                            final int checkPosition;
                            if (expectedSamRecord.getReadNegativeStrandFlag()) {
                                checkPosition = readLength - seqvar.getReadIndex();
                            } else {
                                checkPosition = seqvar.getReadIndex() + i - 1;
                            }
                            compareField("Quality at location specified by seqvar (" + checkPosition + ") doesn't match ",
                                    expectedSamRecord.getBaseQualityString().substring(checkPosition, checkPosition + 1),
                                    actualSamRecord.getBaseQualityString().substring(checkPosition, checkPosition + 1));
                            i++;
                        }
                    }
                }
            }
        }
        readNum++;
        if (!comparisonFailures.isEmpty()) {
            dumpComparison();
            comparisonFailureCount++;
        }
        return comparisonFailures.isEmpty();
    }

    /**
     * Dump the details of expectedSamRecord and actualSamRecord (and gobyAlignment if available). This is
     * called when there are differences between expected and actual to help debug the conversion process.
     * The differences are written to stdout.
     */
    public void dumpComparison() {
        comparisonErrorDump.setLength(0);
        comparisonErrorDump.append("Read Num         : ").append(readNum).append('\n');
        comparisonErrorDump.append("     ERROR(s)    : ").append(ArrayUtils.toString(comparisonFailures)).append('\n');
        if (gobyAlignment != null) {
            comparisonErrorDump.append("     g.index     : ").append(gobyAlignment.getQueryIndex()).append('\n');
            comparisonErrorDump.append("     g.position  : ").append(gobyAlignment.getPosition()).append('\n');
            comparisonErrorDump.append("     g.leftClip  : ").append(gobyAlignment.getSoftClippedBasesLeft()).append('\n');
            comparisonErrorDump.append("     g.rightClip : ").append(gobyAlignment.getSoftClippedBasesRight()).append('\n');
            comparisonErrorDump.append("     g.qAlignLen : ").append(gobyAlignment.getQueryAlignedLength()).append('\n');
            comparisonErrorDump.append("     g.tAlignLen : ").append(gobyAlignment.getTargetAlignedLength()).append('\n');
        }
        comparisonErrorDump.append("     readName (S): ").append(expectedSamRecord.getReadName()).append('\n');
        comparisonErrorDump.append("     readName (D): ").append(actualSamRecord.getReadName()).append('\n');
        comparisonErrorDump.append("     position (S): ").append(expectedSamRecord.getAlignmentStart()).append('\n');
        comparisonErrorDump.append("     position (D): ").append(actualSamRecord.getAlignmentStart()).append('\n');
        comparisonErrorDump.append("     refName  (S): ").append(expectedSamRecord.getReferenceName()).append('\n');
        comparisonErrorDump.append("     refName  (D): ").append(actualSamRecord.getReferenceName()).append('\n');
        comparisonErrorDump.append("     flags    (S): ").append(expectedSamRecord.getFlags()).append('\n');
        comparisonErrorDump.append("     flags    (D): ").append(actualSamRecord.getFlags()).append('\n');
        comparisonErrorDump.append("     mapQual  (S): ").append(expectedSamRecord.getMappingQuality()).append('\n');
        comparisonErrorDump.append("     mapQual  (D): ").append(actualSamRecord.getMappingQuality()).append('\n');
        comparisonErrorDump.append("     negStrand(S): ").append(expectedSamRecord.getReadNegativeStrandFlag()).append('\n');
        comparisonErrorDump.append("     negStrand(D): ").append(actualSamRecord.getReadNegativeStrandFlag()).append('\n');
        comparisonErrorDump.append("     cigar    (S): ").append(expectedSamRecord.getCigar()).append('\n');
        comparisonErrorDump.append("     cigar    (D): ").append(actualSamRecord.getCigar()).append('\n');
        comparisonErrorDump.append("     mdz      (S): ").append(expectedSamRecord.getStringAttribute("MD")).append('\n');
        comparisonErrorDump.append("     mdz      (D): ").append(actualSamRecord.getStringAttribute("MD")).append('\n');
        comparisonErrorDump.append("     mdz-c    (S): ").append(SamHelper.canonicalMdz(expectedSamRecord.getStringAttribute("MD"))).append('\n');
        comparisonErrorDump.append("     mdz-c    (D): ").append(SamHelper.canonicalMdz(actualSamRecord.getStringAttribute("MD"))).append('\n');
        comparisonErrorDump.append("     read     (S): ").append(expectedSamRecord.getReadString()).append('\n');
        comparisonErrorDump.append("     read     (D): ").append(actualSamRecord.getReadString()).append('\n');
        comparisonErrorDump.append("     qual     (S): ").append(expectedSamRecord.getBaseQualityString()).append('\n');
        comparisonErrorDump.append("     qual     (D): ").append(actualSamRecord.getBaseQualityString()).append('\n');
        System.out.println(comparisonErrorDump.toString());
    }

    /**
     * Compare an int field.
     * @param error the error message if the comparison failes
     * @param expected the expected value
     * @param actual the actual value
     */
    private void compareField(final String error, final int expected, final int actual) {
        if (expected != actual) {
            comparisonFailures.add(error);
        }
    }

    /**
     * Compare a string field.
     * @param error the error message if the comparison failes
     * @param expected the expected value
     * @param actual the actual value
     */
    private void compareField(final String error, final String expected, final String actual) {
        if (!expected.equals(actual)) {
            comparisonFailures.add(error);
        }
    }

    /**
     * Compare a boolean field.
     * @param error the error message if the comparison failes
     * @param expected the expected value
     * @param actual the actual value
     */
    private void compareField(final String error, final boolean expected, final boolean actual) {
        if (expected != actual) {
            comparisonFailures.add(error);
        }
    }
}
