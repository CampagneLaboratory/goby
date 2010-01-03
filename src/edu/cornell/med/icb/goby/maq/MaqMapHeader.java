/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.goby.maq;

import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.lang.ArrayUtils;

import java.util.List;

/**
 * Class to store a MAQ MAP header.
 *
 * @author Kevin Dorff
 */
public class MaqMapHeader {
    /**
     * The format.
     */
    private int format;           //int

    // nRef (int) is missing as it is no longer needed since
    // we are using a List[String] for refName

    private IndexedIdentifier referenceNames;
    private DoubleIndexedIdentifier indexToReferenceName;

    /**
     * The tranlsation array. Used when translating from one MaqMapHeader to another.
     * This is generally not used, but classes like FaxMaqMapsRefsMode will use this.
     * The header we are translating TO should be the same size or smaller than the
     * one we are translating from.
     */
    private int[] translation;

    /**
     * The number of reads, may be zero.
     */
    private long numberOfReads;

    /**
     * Constructor.
     */
    public MaqMapHeader() {
        referenceNames = new IndexedIdentifier();
        format = MaqConstants.MAQMAP_FORMAT_NEW;
        translation = null;
    }

    /**
     * Get format.
     *
     * @return format
     */
    public int getFormat() {
        return format;
    }

    /**
     * Set format.
     *
     * @param format format
     */
    public void setFormat(final int format) {
        this.format = format;
    }

    /**
     * Get nRef.
     *
     * @return nRef
     */
    public int getNRef() {
        return referenceNames.size();
    }

    /**
     * Get all refName values.
     *
     * @return all refName values.
     */
    public List<String> getRefName() {
        final int numRefNames = getNRef();
        final List<String> list = new ObjectArrayList<String>(numRefNames);
        for (int i = 0; i < numRefNames; i++) {
            list.add(getRefName(i));
        }
        return list;
    }

    /**
     * Clear the refNames so they can be re-added. DO NOT USE THIS IF YOU DON'T KNOW
     * WHAT ARE YOU DOING!! WHen you add the references back in (using addRefName()
     * take care to add them in the right order).
     */
    public void clearRefNames() {
        this.referenceNames = new IndexedIdentifier();
        this.indexToReferenceName = null;
    }

    /**
     * Add a refName value. If newRefName already exists in refName, this will NOT add it again.
     *
     * @param newRefName the new refName to add
     * @return the position (0-based) within refName where newRefName is stored
     */
    public synchronized int addRefName(final String newRefName) {
        return referenceNames.registerIdentifier(new MutableString(newRefName));
    }

    /**
     * Get the index of a specific refName.
     *
     * @param existingRefName the refName to get the index for
     * @return the refName index position or -1 if not found
     */
    public synchronized int getRefNameIndex(final String existingRefName) {
        return referenceNames.getInt(new MutableString(existingRefName));
    }

    /**
     * Get one of the refName values.
     *
     * @param i the refName index to get
     * @return one of the refName values
     */
    public String getRefName(final int i) {
        if (indexToReferenceName == null || indexToReferenceName.size() != referenceNames.size()) {
            referencesDefined();
        }

        final MutableString id = this.indexToReferenceName.getId(i);
        return id == null ? null : id.toString();
    }

    /**
       * Get one of the refName values.
       *
       * @param i the refName index to get
       * @return one of the refName values
       */
    public MutableString getRefNameMutable(final int i) {
        if (indexToReferenceName == null || indexToReferenceName.size() != referenceNames.size()) {
            referencesDefined();
        }

        return this.indexToReferenceName.getId(i);
    }

    /**
     * Get numberOfReads (may be zero).
     *
     * @return numberOfReads
     */
    public long getNumberOfReads() {
        return numberOfReads;
    }

    /**
     * Set numberOfReads. This is often just 0 and cannot be relied upon.
     *
     * @param numberOfReads numberOfReads
     */
    public void setNumberOfReads(final long numberOfReads) {
        this.numberOfReads = numberOfReads;
    }



    /**
     * Copy this object to a new object.
     *
     * @return a copy of this object
     */
    public MaqMapHeader copy() {
        final MaqMapHeader result = new MaqMapHeader();
        result.format = format;
        result.referenceNames = new IndexedIdentifier();
        for (final MutableString name : this.referenceNames.keySet()) {
            result.referenceNames.put(name, this.referenceNames.get(name));
        }
        result.referencesDefined();
        result.numberOfReads = numberOfReads;
        return result;
    }

    public void referencesDefined() {
        this.indexToReferenceName = new DoubleIndexedIdentifier(referenceNames);
    }

    /**
     * Merge header information from another MaqMapHeader into this one.
     *
     * @param other the other MaqMapHeader
     */
    public void merge(final MaqMapHeader other) {
        for (final String currentRefName : other.getRefName()) {
            addRefName(currentRefName);
        }
        numberOfReads += other.getNumberOfReads();
    }

    /**
     * Human readable information about this object.
     *
     * @return Human readable information about this object
     */
    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append("format=").append(format).append(", ");
        sb.append("n_ref=").append(getNRef()).append(", ");
        sb.append("numberOfReads=").append(numberOfReads).append('\n');
        for (int i = 0; i < getNRef(); i++) {
            sb.append("refName[").append(i).append("]=");
            sb.append(ArrayUtils.toString(getRefName(i))).append("\n");
        }
        return sb.toString();
    }

    /**
     * Given the index in the old header, return the index in the new header.
     * @param index the index in the old header
     * @return the index in the new header
     */
    public int getTranslationForIndex(final int index) {
        if (translation != null) {
            return translation[index];
        } else {
            return -1;
        }
    }

    /**
     * Set the translation array. Generally not used, but can help with translating
     * refName positions from a different header. The other header should be the
     * same size or larger, translation[] will be the size of the OTHER list of
     * headers and contain the indexes into THIS headers ref names.
     * @param translation the translation array.
     */
    public void setTranslation(final int[] translation) {
        this.translation = translation;
    }
}
