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

import java.util.List;

/**
 * Class for comparing two SAMRecords. Used with both tests and
 * with SAMComparisonMode.
 */
public class SamComparison {

    public boolean mappedQualitiesPreserved;
    public boolean softClipsPreserved;
    public boolean checkMate = true;

    public int readNum;
    public SAMRecord expectedSamRecord;
    public SAMRecord actualSamRecord;
    public Alignments.AlignmentEntry gobyAlignment;
    public int comparisonFailureCount;

    private final MutableString comparisonErrorDump = new MutableString();
    private boolean comparisonFailure;

    public void reset() {
        readNum = 0;
        comparisonFailureCount = 0;
    }

    /**
     * Compare expectedSamRecord vs actualSamRecord.
     * @return true if the two records are found to be the same
     */
    public boolean compareSamRecords() {
        comparisonFailure = false;
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
        compareField("MD:Z doesn't match",
                SamHelper.canonicalMdz(expectedSamRecord.getStringAttribute("MD")),
                SamHelper.canonicalMdz(actualSamRecord.getStringAttribute("MD")));

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
        if (comparisonFailure) {
            comparisonFailureCount++;
        }
        return !comparisonFailure;
    }

    public void dumpComparison(final String error) {
        comparisonErrorDump.setLength(0);
        comparisonErrorDump.append("Read Num         : ").append(readNum).append('\n');
        comparisonErrorDump.append("     ERROR       : ").append(error).append('\n');
        if (gobyAlignment != null) {
            comparisonErrorDump.append("     g.position  : ").append(gobyAlignment.getPosition()).append('\n');
            comparisonErrorDump.append("     g.leftClip  : ").append(gobyAlignment.getSoftClippedBasesLeft()).append('\n');
            comparisonErrorDump.append("     g.rightClip : ").append(gobyAlignment.getSoftClippedBasesRight()).append('\n');
            comparisonErrorDump.append("     g.qAlignLen : ").append(gobyAlignment.getQueryAlignedLength()).append('\n');
            comparisonErrorDump.append("     g.tAlignLen : ").append(gobyAlignment.getTargetAlignedLength()).append('\n');
        }
        comparisonErrorDump.append("     position (S): ").append(expectedSamRecord.getAlignmentStart()).append('\n');
        comparisonErrorDump.append("     position (S): ").append(expectedSamRecord.getAlignmentStart()).append('\n');
        comparisonErrorDump.append("     refName  (S): ").append(expectedSamRecord.getReferenceName()).append('\n');
        comparisonErrorDump.append("     refName  (D): ").append(actualSamRecord.getReferenceName()).append('\n');
        comparisonErrorDump.append("     readName (S): ").append(expectedSamRecord.getReadName()).append('\n');
        comparisonErrorDump.append("     readName (D): ").append(actualSamRecord.getReadName()).append('\n');
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


    private void compareField(final String error, final int expected, final int actual) {
        if (!comparisonFailure && expected != actual) {
            comparisonFailure = true;
            dumpComparison(error);
        }
    }

    private void compareField(final String error, final String expected, final String actual) {
        if (!comparisonFailure && !expected.equals(actual)) {
            comparisonFailure = true;
            dumpComparison(error);
        }
    }

    private void compareField(final String error, final boolean expected, final boolean actual) {
        if (!comparisonFailure && expected != actual) {
            comparisonFailure = true;
            dumpComparison(error);
        }
    }
}
