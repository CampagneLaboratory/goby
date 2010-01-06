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

import org.apache.commons.configuration.BaseConfiguration;
import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.PropertiesConfiguration;
import org.apache.commons.configuration.SystemConfiguration;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.Iterator;

/**
 * Helper class to load project-wide properties.
 *
 * @see GobyPropertyKeys
 * @author Fabien Campagne
 *         Date: Jul 26, 2009
 *         Time: 4:16:33 PM
 */
public class ConfigHelper {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(ConfigHelper.class);

    /**
     * Configuration properties used by goby.
     */
    private final CompositeConfiguration configuration;

    /**
     * List of files to configure goby with.  Note that the files will by loaded in order
     * and the first one found will be used.
     */
    private static final String[] DEFAULT_CONFIG_FILE_LOCATIONS = {
            "goby.properties",
            "config/goby.properties"
    };

    /**
     * Singleton instance of this class.
     */
    private static final ConfigHelper INSTANCE = new ConfigHelper(DEFAULT_CONFIG_FILE_LOCATIONS);

    /**
     * Load the Goby configuration.
     * @param defaultConfigFileLocations locations for configurations to check if one
     * was not explicitly defined in a system property.
     */
    private ConfigHelper(final String... defaultConfigFileLocations) {
        super();
        configuration = new CompositeConfiguration();

        // by loading the system configuration first, any values set on the java command line
        // with "-Dproperty.name=property.value" will take precedence
        configuration.addConfiguration(new SystemConfiguration());

        // check to see if the user specified a configuration in the style of log4j
        if (configuration.containsKey("goby.configuration")) {
            final String configurationURLString = configuration.getString("goby.configuration");
            try {
                if (LOG.isDebugEnabled()) {
                    LOG.debug("Attempting to load " + configurationURLString);
                }
                final URL configurationURL = new URL(configurationURLString);
                configuration.addConfiguration(new PropertiesConfiguration(configurationURL));
                LOG.info("Goby configured with " + configurationURL);
            } catch (MalformedURLException e) {
                LOG.error("Invalid Goby configuration", e);
            } catch (ConfigurationException e) {
                LOG.error("Could not configure Goby from " + configurationURLString, e);
            }
        } else {
            // no configuration file location specified so check the default locations
            for (final String configFile : defaultConfigFileLocations) {
                try {
                    if (LOG.isDebugEnabled()) {
                        LOG.debug("Attempting to load " + configFile);
                    }
                    configuration.addConfiguration(new PropertiesConfiguration(configFile));
                } catch (ConfigurationException e) {
                    continue;
                }

                // if we got here the file was loaded and we don't search any further
                LOG.info("Goby configured with " + configFile);
                break;
            }
        }

        // load "default" configurations for any properties not found elsewhere
        // it's important that this added LAST so the user can override any settings
        final Configuration defaultConfiguration = getDefaultConfiguration();
        configuration.addConfiguration(defaultConfiguration);

        if (LOG.isDebugEnabled()) {
            LOG.debug("Goby configuration: ");
            final Iterator keys = configuration.getKeys();
            while (keys.hasNext()) {
                final String key = (String) keys.next();
                LOG.debug(key + " = " + configuration.getString(key));
            }
        }
    }

    /**
     * Any default configuration item values if required should defined here.
     * @return The default configuration to use.
     */
    private Configuration getDefaultConfiguration() {
        final Configuration defaultConfiguration = new BaseConfiguration();
        defaultConfiguration.addProperty(GobyPropertyKeys.EXECUTABLE_PATH_LASTAG, ".");
        defaultConfiguration.addProperty(GobyPropertyKeys.EXECUTABLE_PATH_BWA, ".");
        defaultConfiguration.addProperty(GobyPropertyKeys.DATABASE_DIRECTORY, ".");
        defaultConfiguration.addProperty(GobyPropertyKeys.WORK_DIRECTORY, ".");
        return defaultConfiguration;
    }

    /**
     * The configuration items for Goby.
     * @return The configuration needed for goby operations.
     */
    public static Configuration getConfiguration() {
        return INSTANCE.configuration;
    }
}
