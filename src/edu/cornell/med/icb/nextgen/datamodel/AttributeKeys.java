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

import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;

/**
 * Keeps default keys together.
 *
 * @author Fabien Campagne
 *         Date: May 21, 2009
 *         Time: 12:04:07 PM
 */
public class AttributeKeys {
    private ObjectSet<String> suggestedKeys;
    private ObjectSet<String> requiredKeys;

    public AttributeKeys() {
        requiredKeys = new ObjectOpenHashSet<String>();
        suggestedKeys = new ObjectOpenHashSet<String>();
    }

    public ObjectSet<String> getSuggestedKeys() {
        return suggestedKeys;
    }

    public void addSuggestedKey(final String name) {
        suggestedKeys.add(name);
    }

    public ObjectSet<String> getRequiredKeys() {
        return requiredKeys;
    }

    public void addRequiredKey(final String name) {
        requiredKeys.add(name);
    }
}
