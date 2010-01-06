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

import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Iterator;

/**
 * Helper class to load project-wide properties.
 *
 * @author Fabien Campagne
 *         Date: Jul 26, 2009
 *         Time: 4:16:33 PM
 */
public class GobyConfiguration {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(GobyConfiguration.class);

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
    private static final GobyConfiguration INSTANCE =
            new GobyConfiguration(DEFAULT_CONFIG_FILE_LOCATIONS);

    /**
     * Path to the directory that contains the lastag executable.
     */
    public static final String EXECUTABLE_PATH_LASTAG = "executables.path.lastag";

    /**
     * Path to the directory that contains the BWA executable.
     */
    public static final String EXECUTABLE_PATH_BWA = "executables.path.bwa";

    /**
     * Path to the work directory. This should be a large scratch location, where results of
     * intermediate calculations will be stored.
     */
    public static final String WORK_DIRECTORY = "work.directory";

    /**
     * Path to the database directory. This is the location where the indexed database files may
     * reside.
     */
    public static final String DATABASE_DIRECTORY = "database.directory";

    /**
     * Load the Goby configuration.
     * @param defaultConfigFileLocations locations for configurations to check if one
     * was not explicitly defined in a system property.
     */
    private GobyConfiguration(final String... defaultConfigFileLocations) {
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
                LOG.info("Goby configured with " + new File(configFile).getAbsolutePath());
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
        defaultConfiguration.addProperty(EXECUTABLE_PATH_LASTAG, "");
        defaultConfiguration.addProperty(EXECUTABLE_PATH_BWA, "");
        defaultConfiguration.addProperty(DATABASE_DIRECTORY, ".");
        defaultConfiguration.addProperty(WORK_DIRECTORY, ".");
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
