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

import edu.cornell.med.icb.io.ResourceFinder;
import groovy.lang.Closure;
import groovy.util.ConfigObject;
import groovy.util.ConfigSlurper;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

/**
 * This is a class that can handle either normal Java properties files or more
 * sophisticated Groovy Config files.
 *
 * @author Kevin Dorff
 */
public class GroovyProperties {

    /**
     * The logger.
     */
    private static final Log LOG = LogFactory.getLog(GroovyProperties.class);

    /**
     * The actual properties are stored here.
     */
    private final Map<String, Object> properties = new HashMap<String, Object>();

    /**
     * Private constructor. Use ConfigHelper to load the confiruation.
     */
    private GroovyProperties() {
    }

    /**
     * Load properties for the specified platforms from the config and/or property files specified.
     * This will try to load each of each platform within configFileLocations. If multiple
     * platforms and multiple config files are specified, the order it will read is
     * - platform1 / file1
     * - platform2 / file1
     * - platform1 / file2
     * - platform2 / file2
     * Any duplication configuration will be merged with the later read configuration
     * taking precidence. This will seamlessly handly both Groovy Config files and
     * Java .properties files. If No configFileLocations are specified (if no files
     * are read) this will read Java's System.Properties.
     *
     * @param platforms           the platform or platforms to load (in order). If more than one
     *                            platform is specified, they will be merged.
     * @param configFilePaths This will first search the class path and then these paths for the
     * specified config files.
     * @param configFileNames the config file names to look for within the classpath and within
     * configFilePaths
     * @return the GroovyProperties
     */
    public static GroovyProperties load(
            final String[] platforms,
            final String[] configFilePaths,
            final String... configFileNames) {

        final GroovyProperties groovyProperties = new GroovyProperties();
        final ResourceFinder resourceFinder = new ResourceFinder(configFilePaths);

        final List<ConfigObject> configObjects = new LinkedList<ConfigObject>();
        if (configFileNames != null) {
            for (final String configFileName : configFileNames) {
                if (configFileName == null) {
                    continue;
                }
                final URL resourceUrl = resourceFinder.findResource(configFileName);
                if (resourceUrl != null) {
                    for (final String platform : platforms) {
                        final ConfigObject configObject = loadConfig(platform, resourceUrl);
                        if (configObject != null) {
                            configObjects.add(configObject);
                        }
                    }
                }
            }
        }
        if (configObjects.isEmpty()) {
            // Load from System.properties since no other configs were loaded
            configObjects.add(loadConfig(null, null));
        }

        // Load the config into our properties map and merge as we go
        for (final ConfigObject currentConfig : configObjects) {
            final Map origProperties = currentConfig.flatten();
            for (final Object key : origProperties.keySet()) {
                groovyProperties.properties.put(key.toString().intern(), origProperties.get(key));
            }
        }
        return groovyProperties;
    }

    /**
     * Load a siingle configuration for a single platform, file.
     *
     * @param platform    the platform to load the configuration file
     * @param configToUse the configuration file to load
     * @return the loaded configuration
     */
    private static ConfigObject loadConfig(final String platform, final URL configToUse) {
        // Support normal properties files, too, here
        String configUrlString = "";
        if (configToUse != null) {
            configUrlString = configToUse.toString();
            LOG.info("Loading configuration from URL " + configUrlString);
        }
        if (configToUse == null || configUrlString.endsWith(".properties")) {
            // system properties or normal properties file
            Properties p;
            if (configToUse != null) {
                p = new Properties();
                InputStream propsInputStream = null;
                try {
                    propsInputStream = configToUse.openStream();
                    LOG.info("Loading java properties from " + configToUse);
                    p.load(propsInputStream);
                } catch (IOException e) {
                    LOG.error("Couldn't load properties from " + configToUse, e);
                    LOG.info("Loading System.properties");
                    p = System.getProperties();
                } finally {
                    IOUtils.closeQuietly(propsInputStream);
                }
            } else {
                LOG.info("Loading System.properties");
                p = System.getProperties();
            }
            final ConfigObject configObject = new ConfigObject();
            for (final String key : p.stringPropertyNames()) {
                final String internedKey = key.intern();
                configObject.put(internedKey, p.getProperty(internedKey));
            }
            return configObject;
        } else {
            return new ConfigSlurper(platform).parse(configToUse);
        }
    }

    /**
     * Translate key values. Good for removing parts of the key, changing "_" to ".", etc.
     *
     * @param fromValue the old value to look for
     * @param toValue   the value to replace the old value with
     */
    public void translateKeys(final String fromValue, final String toValue) {
        final List<String> keys = new LinkedList<String>();
        keys.addAll(properties.keySet());
        for (String oldKey : keys) {
            oldKey = oldKey.intern();
            if (oldKey.contains(fromValue)) {
                final String newKey = StringUtils.replace(oldKey, fromValue, toValue);
                final Object value = properties.get(oldKey);
                properties.remove(oldKey);
                properties.put(newKey.intern(), value);
            }
        }
    }

    /**
     * Get the value for the specific key. If the value isn't found, null will be returned.
     *
     * @param key the key to look for
     * @return the value for that key
     */
    public String get(final String key) {
        return get(key, null);
    }

    /**
     * Get the value for the specific key. If the key isn't found, defaultValue will be returned.
     * IF key is assocaited with a closure, an empty string will be passed as the parameter
     * to that closure.
     *
     * @param key          the key to look for
     * @param defaultValue the default value to return if key isn't assocated with a value
     * @return the value for that key
     */
    public String get(final String key, final String defaultValue) {
        final Object value = properties.get(key.intern());
        if (value == null) {
            return defaultValue;
        }
        if (value instanceof Closure) {
            return getWithParameter(key, "", defaultValue);
        } else {
            return value.toString();
        }
    }


        /**
         * Get the value for the specific key with the idea that the key is probably associated
         * with a closure. If the key isn't found, null will be returned.
         * If key is associated with a closure, closureParameter will be passed as a parameter
         * to the closure, if key isn't associated with a closure the associated value will be returned.
         *
         * @param key              the key to look for
         * @param closureParameter the paramter to be passed to the closure (if key is associated
         *                         with a closure).
         * @return the value for that key
         */

    public String getWithParameter(final String key, final String closureParameter) {
        return getWithParameter(key, closureParameter, null);
    }

    /**
     * Get the value for the specific key with the idea that the key is probably assocaited
     * with a closure. If the key isn't found, defaultValue will be returned.
     * If key is associated with a closure, closureParameter will be passed as a parameter
     * to the closure, if key isn't associated with a closure the associated value will be returned.
     *
     * @param key              the key to look for
     * @param closureParameter the paramter to be passed to the closure (if key is associated
     *                         with a closure).
     * @param defaultValue     the default value to return if key isn't assocated with a value
     * @return the value for that key
     */
    public String getWithParameter(
            final String key, final String closureParameter, final String defaultValue) {
        final Object value = properties.get(key.intern());
        if (value == null) {
            return defaultValue;
        }
        if (value instanceof Closure) {
            final Closure closure = (Closure) value;
            return closure.call(closureParameter).toString();
        } else {
            return value.toString();
        }
    }

    /**
     * Return the set of keys, the names of the properties stored in this object.
     *
     * @return the set of keys
     */
    public Set<String> keySet() {
        final Set<String> copy = new LinkedHashSet<String>();
        copy.addAll(properties.keySet());
        return copy;
    }

    /**
     * Display the properties as a String.
     *
     * @return properties as string
     */
    @Override
    public String toString() {
        return ArrayUtils.toString(properties);
    }

    public String assertGet(final String key) {
        final String value = get(key);
        assert value != null : "Property " + key + " must be defined in goby config.";
        return value;
    }


    public String assertGetWithParameter(final String key, final String parameter) {
        final String value = getWithParameter(key, parameter);
        assert value != null : "Property " + key + " must be defined in goby config.";
        return value;
    }

    public String assertGetWithParameter(final String key, final String parameter, final String defaultValue) {
        final String value = getWithParameter(key, parameter, defaultValue);
        assert value != null : "Property " + key + " must be defined in goby config.";
        return value;
    }
}
