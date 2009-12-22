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

package edu.cornell.med.icb.util;

import org.junit.Test;
import static org.junit.Assert.assertEquals;

import java.io.IOException;

/**
 * Test loading groovy config and java properties files into GroovyProperties objects.
 * @author Kevin Dorff
 */
public class TestGroovyProperties {
    @Test
    public void testLoadProperties() throws IOException {
        GroovyProperties props = GroovyProperties.load(new String[] {"dubble"}, new String[] {"test-data"}, "TestConfig.groovy", null);
        System.out.println(props);
        props.translateKeys("settings.", "");
        assertEquals("1", props.get("one"));
        assertEquals("8", props.get("four"));
        assertEquals("4", props.get("two"));
        assertEquals(null, props.get("five"));
        assertEquals("five", props.get("five", "five"));
        assertEquals("two-5", props.getWithParameter("any", "5"));
        assertEquals("none", props.getWithParameter("anyx", "5", "none"));
    }

    @Test
    public void testLoadMultipePropertiesMultipleLocations() throws IOException {
        GroovyProperties props = GroovyProperties.load(new String[] {"dubble", "quadruple_four"}, new String[] {"test-data"}, "blah", "TestConfig.groovy");
        assertEquals("1", props.get("settings.one"));
        assertEquals("sixteen", props.get("settings.four"));
        assertEquals("4", props.get("settings.two"));
        assertEquals(null, props.get("settings.five"));
        assertEquals("two-5", props.getWithParameter("settings.any", "5"));
        assertEquals("none", props.getWithParameter("settings.anyx", "5", "none"));
    }

    @Test
    public void testNormalProperties() throws IOException {
        // Two valid files are presented, it will use the first one it files which happens to be a properties file
        GroovyProperties props = GroovyProperties.load(new String[] {"dubble"}, new String[] {"test-data"}, "TestConfig.groovy", "test.properties");
        props.translateKeys("settings.", ""); // Won't do anything
        assertEquals("1", props.get("one"));
        assertEquals("two", props.get("two"));
        assertEquals("three", props.get("three"));
        assertEquals("8", props.get("four"));
        assertEquals(null, props.get("five"));
        assertEquals("5", props.get("five", "5"));
    }

}
