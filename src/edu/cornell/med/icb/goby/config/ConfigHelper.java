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

package edu.cornell.med.icb.goby.config;

import edu.cornell.med.icb.goby.util.GroovyProperties;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.LinkedList;
import java.util.List;

/**
 * Helper class to load project-wide properties.
 *
 * @author Fabien Campagne
 *         Date: Jul 26, 2009
 *         Time: 4:16:33 PM
 */
public class ConfigHelper {
    private static final Log LOG = LogFactory.getLog(ConfigHelper.class);
    private static GroovyProperties singleton;

    private ConfigHelper() {
        super();
    }

    public static GroovyProperties loadConfiguration() {
        return loadConfiguration(
                new String[] {"config/", "../config/", "../../config/", "./"},
                "GlobalConfig.groovy",
                "LocalConfig.groovy");
    }

    /**
     * Load properties or return the JVM-wide singleton if already loaded.
     * This version gives full control over the config directories and config
     * filenames for global and local.
     * @return
     */
    public static GroovyProperties loadConfiguration(
            final String[] configDirs,
            final String globalConfigFilename,
            final String localConfigFilename) {

        if (singleton != null) {
            return singleton;
        }

        final GroovyProperties properties = GroovyProperties.load(
                obtainPlatforms(), configDirs, globalConfigFilename, localConfigFilename);
        properties.translateKeys("settings.", "");
        properties.translateKeys("_", ".");
        singleton = properties;

        if (LOG.isDebugEnabled()) {
            LOG.debug("Loaded Config=" + properties.toString());
        }

        return properties;
    }

    /**
     * Normalize an array of dirnames. Makes the separators unix and makes sure
     * each directory name ends in "/". dirs is not changed in place, a new array
     * is made and returned (unles dirs is empty).
     * @param dirs the directories
     * @return the directories normalized.
     */
    static String[] normalizeDirs(final String[] dirs) {
        if (dirs == null || dirs.length == 0) {
            return dirs;
        }
        final String[] result = new String[dirs.length];
        for (int i = 0; i < dirs.length; i++) {
            final String dir = FilenameUtils.separatorsToUnix(dirs[i]);
            if (dir.endsWith("/")) {
                result[i] = dir;
            } else {
                result[i] = dir + "/";
            }
        }
        return result;
    }

    /**
     * Expands file with each of the directory names. If dirs is empty or null
     * this will return just file.
     * @param dirs assumed to be a list of dirs with unix separators ending in "/"
     * (can call normalizeDirs to fix them before calling this).
     * @param file the file to expand
     * @return the expanded filenames including the directories
     */
    static String[] expandDirs(final String[] dirs, final String file) {
        if (file == null) {
            return null;
        }
        if (dirs == null || dirs.length == 0) {
            return new String[] {file};
        }
        final String[] result = new String[dirs.length];
        for (int i = 0; i < dirs.length; i++) {
            result[i] = dirs[i] + file;
        }
        return result;
    }

    /**
     * Determine the array of platforms to load (for Config.groovy).
     * (one or more of amazon, cornell, windows).
     *
     * @return the array of platforms
     */
    private static String[] obtainPlatforms() {
        final List<String> configList = new LinkedList<String>();
        String addressLastPiece;
        try {
            final InetAddress localHost = InetAddress.getLocalHost();
            final String[] parts = localHost.getCanonicalHostName().split("\\.");
            addressLastPiece = parts[parts.length - 1];
        } catch (UnknownHostException e) {
            LOG.error("Couldn't not obtain InetAddress.getLocalHost(). Assuming .edu address.", e);
            addressLastPiece = "edu";
        }

        if ("edu".equals(addressLastPiece)) {
            configList.add("cornell");
        } else if ("cluster".equals(addressLastPiece)) {
            configList.add("PBS");
        } else {
            configList.add("amazon");
        }
        return configList.toArray(new String[configList.size()]);
    }
}
