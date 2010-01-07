/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

package edu.cornell.med.icb.goby.util;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Ability to write data that goes to the OutputStream to the log at
 * a specified log level.
 *
 * @author Kevin Dorff
 */
public class LoggingOutputStream extends OutputStream {

    /** Store data until end of line. */
    private final ByteArrayOutputStream buffer = new ByteArrayOutputStream();

    /** Start of line prefixes of output. */
    private byte[] prefix;

    /** Start of line prefixes of output. */
    private final String prefixStr;

    /** The log level to write logs at. */
    private final Level logLevel;

    /** The logger to log to. */
    private final Logger log;

    /**
     * Construct.
     * @param logAsClass the class to log as
     * @param logLevel the log level to log at
     * @param prefixStr the optional prefix that will be put at the front of each log line.
     */
    public LoggingOutputStream(
            final Class logAsClass, final Level logLevel, final String prefixStr) {
        log = Logger.getLogger(logAsClass);
        if (StringUtils.isBlank(prefixStr)) {
            this.prefix = null;
            this.prefixStr = null;
        } else {
            this.prefixStr = prefixStr;
            this.prefix = prefixStr.getBytes();
        }
        this.logLevel = logLevel;
    }

    /**
     * Replace __OUTPUT_TAG__ in prefix with this output tag.
     * @param outputTag the tag to put in the prefix (if the prefix contains "__OUTPUT_TAG__")
     */
    final void setOutputTag(final String outputTag) {
        if (prefixStr == null) {
            return;
        }
        final String newPrefixStr = StringUtils.replace(prefixStr, "__OUTPUT_TAG__", outputTag);
        this.prefix = newPrefixStr.getBytes();
    }

    /** Last character received. */
    int lastb = 0;

    /**
     * Write a character to the output stream. At the end of the line this will log the output.
     * @param b the character written to the output stream.
     */
    public void write(int b) {
        if (b == 13) {
            // Trace CR as LF
            b = 10;
        }
        if (b == 10) {
            if (lastb == 10) {
                // Ignore blank lines (and if we get CR/LF treat it as one character instead of 2)
            } else {
                log.log(logLevel, buffer.toString());
                buffer.reset();
                try {
                    if (prefix != null) {
                        buffer.write(prefix);
                    }
                } catch (IOException e) {
                    log.error("Cannot prepend prefix to buffer");
                }
            }
        } else {
            buffer.write((char) b);
        }
        lastb = b;
    }
}
