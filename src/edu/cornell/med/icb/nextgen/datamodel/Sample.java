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

package edu.cornell.med.icb.nextgen.datamodel;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;

import java.util.Collections;

/**
 * Represents a next-gen sample. The attributes platform and organism must be defined for each sample.
 *
 * @author Fabien Campagne
 *         Date: May 21, 2009
 *         Time: 11:47:56 AM
 */
public class Sample extends Tagged {
    private static AttributeKeys attributeKeys;
    static {
        attributeKeys = new AttributeKeys();
        attributeKeys.addRequiredKey("platform");
        attributeKeys.addRequiredKey("organism");

        attributeKeys.addSuggestedKey("tissue");
        attributeKeys.addSuggestedKey("protocol");
    }

    /**
     * Short name for this sample.
     */
    String name;

    /**
     * Sample description, unstructured text.
     */
    String description;

    /**
     * Where the file uploaded by the user is stored.
     */
    StoredFile upload;

    /**
     * Where the files converted to compact-reads format are stored. Null when the upload has not yet been converted.
     */
    StoredFile[] compactReads;

    /**
     * The number of reads in the sample.
     */
    int numberOfReads;

    /**
     * The minimum read length in this sample.
     */
    int minimumLength;

    /**
     * The maximum read length in this sample.
     */
    int maximumLength;

    /**
     * Attributes describe this sample in specific ways. They can be used to organize samples into
     * groups.
     */
    private Attributes attributes;

    public Sample() {
        attributes = new Attributes();
    }

    public static AttributeKeys getAttributeKeys() {
        return attributeKeys;
    }

    public Attributes getAttributes() {
        return attributes;
    }

    public StoredFile[] getCompactReads() {
        return compactReads;
    }

    public void setCompactReads(final StoredFile[] compactReads) {
        this.compactReads = compactReads;
    }

    /**
     * Add a compact reads stored file.
     * @param compactReadsToAdd
     */
    public void addCompactReads(final StoredFile... compactReadsToAdd) {
        final ObjectList<StoredFile> result = new ObjectArrayList<StoredFile>();
        if (this.compactReads != null) {
            Collections.addAll(result, this.compactReads);
        }
        for (final StoredFile file : compactReadsToAdd) {
            result.add(file);
        }
        this.compactReads = result.toArray(new StoredFile[result.size()]);
    }

    public int getMinimumLength() {
        return minimumLength;
    }

    public void setMinimumLength(final int minimumLength) {
        this.minimumLength = minimumLength;
    }

    public int getMaximumLength() {
        return maximumLength;
    }

    public void setMaximumLength(final int maximumLength) {
        this.maximumLength = maximumLength;
    }


    public StoredFile getUpload() {
        return upload;
    }

    public void setUpload(final StoredFile upload) {
        this.upload = upload;
    }

    public int getNumberOfReads() {
        return numberOfReads;
    }

    public void setNumberOfReads(final int numberOfReads) {
        this.numberOfReads = numberOfReads;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(final String description) {
        this.description = description;
    }

    public void setName(final String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }
}
