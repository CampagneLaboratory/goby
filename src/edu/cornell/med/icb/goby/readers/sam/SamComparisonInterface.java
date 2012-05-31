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
import net.sf.samtools.SAMRecord;

/**
 * Interface for comparing SAM records.
 */
public interface SamComparisonInterface {
    public void reset();
    public int compare(final SAMRecord source, final SAMRecord dest, final Alignments.AlignmentEntry gobyDest);
    public int finished();

    /**
     * Get if it is assumed that the compact file created from the BAM/SAM
     * file preserved mapped qualities.
     * @return if it is assumed ...
     */
    public boolean isMappedQualitiesPreserved();

    /**
     * Set if it is assumed that the compact file created from the BAM/SAM
     * file preserved mapped qualities.
     * @param mappedQualitiesPreserved if it is assumed...
     */
    public void setMappedQualitiesPreserved(final boolean mappedQualitiesPreserved);

    /**
     * Get if it is assumed that the compact file created from the BAM/SAM
     * file preserved soft clips.
     * @return if it is assumed ...
     */
    public boolean isSoftClipsPreserved();

    /**
     * Set if it is assumed that the compact file created from the BAM/SAM
     * file preserved soft clips.
     * @param softClipsPreserved if it is assumed ...
     */
    public void setSoftClipsPreserved(final boolean softClipsPreserved);

    /**
     * Get if the details about mate reads will be checked.
     * If the source SAM/BAM file is a complete file you can set this to true,
     * if you are using an incomplete source SAM/BAM file, this should be
     * set to false. Default is false.
     * @return if mates will be checked
     */
    public boolean isCheckMate();

    /**
     * Set if the details about mate reads will be checked.
     * If the source SAM/BAM file is a complete file you can set this to true,
     * if you are using an incomplete source SAM/BAM file, this should be
     * set to false. Default is false.
     * @param checkMate if mates will be checked
     */
    public void setCheckMate(final boolean checkMate);

    /**
     * Get if canonical MD:Z comparisons will be made.
     * When true, the source and destination MD:Z values will be passed through an algorithm
     * to make them canonical (place 0's in places where 0's should exist but might not).
     * By default this is enabled.
     * @return if ...
     */
    public boolean isCanonicalMdzForComparison();

    /**
     * Set if canonical MD:Z comparisons will be made.
     * When true, the source and destination MD:Z values will be passed through an algorithm
     * to make them canonical (place 0's in places where 0's should exist but might not).
     * By default this is enabled.
     * @param canonicalMdzForComparison if ...
     */
    public void setCanonicalMdzForComparison(final boolean canonicalMdzForComparison);

    /**
     * Get if read names were preserved and thus we can use them to aid with comparison.
     * @return if...
     */
    public boolean isReadNamesPreserved();

    /**
     * Set if read names were preserved and thus we can use them to aid with comparison.
     * @param readNamesPreserved if...
     */
    public void setReadNamesPreserved(final boolean readNamesPreserved);

    /**
     * Return how many reads have been compared since reset() was last called.
     * @return how many...
     */
    public int getReadNum();

    /**
     * Return how many comparison failures have been found since reset() was last called.
     * @return how many...
     */
    public int getComparisonFailureCount();
}
