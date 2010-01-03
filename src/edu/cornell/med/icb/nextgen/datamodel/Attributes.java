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

import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;

/**
 *
 * @author Fabien Campagne
 *         Date: May 21, 2009
 *         Time: 12:08:43 PM
 */
class Attributes {
    /**
     * Sample attributes. Used to annotate a sample with its characteristics.
     * key=tisssue, value=brain indicates that the sample was obtained from brain.
     */
    private Object2ObjectMap<String, String> attributes;

    public Attributes() {
        attributes = new Object2ObjectOpenHashMap<String, String>();
    }

    /**
     * Returns true if this sample has this attribute.
     *
     * @param key attribute key.
     * @return True or False.
     */
    public boolean hasAttribute(final String key) {
        return attributes.containsKey(key);

    }

    public String getAttribute(final String key) {
        return attributes.get(key);
    }

    public void setAttribute(final String key, final String value) {
        attributes.put(key, value);
    }
}
