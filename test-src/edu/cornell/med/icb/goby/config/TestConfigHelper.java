/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.config;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertNull;
import org.junit.Test;

/**
 * Describe class here.
 *
 * @author Kevin Dorff
 */
public class TestConfigHelper {
    private final String[] emptyDirsArray = new String[0];
    private final String[] unNormalized = {
            "\\this\\dir\\path", "another/dir\\path", "yet\\another/path\\",
            "onemore/", "finalpath" };
    private final String[] normalized = {
            "/this/dir/path/", "another/dir/path/", "yet/another/path/",
            "onemore/", "finalpath/" };
    private final String[] normalizedExpanded = {
            "/this/dir/path/cat", "another/dir/path/cat", "yet/another/path/cat",
            "onemore/cat", "finalpath/cat" };

    /**
     * Test than directory normalization expansion works as expected.
     */
    @Test
    public void testNormalize() {
        assertNull(ConfigHelper.normalizeDirs(null));
        assertArrayEquals(emptyDirsArray, ConfigHelper.normalizeDirs(emptyDirsArray));
        assertArrayEquals(normalized, ConfigHelper.normalizeDirs(unNormalized));
    }

    /**
     * Test than expansion works as expected.
     */
    @Test
    public void testExpand() {
        assertNull(null, ConfigHelper.expandDirs(null, null));
        assertArrayEquals(new String[] {"cat"}, ConfigHelper.expandDirs(null, "cat"));
        assertArrayEquals(new String[] {"cat"}, ConfigHelper.expandDirs(emptyDirsArray, "cat"));
        assertArrayEquals(normalizedExpanded, ConfigHelper.expandDirs(normalized, "cat"));
    }
}
