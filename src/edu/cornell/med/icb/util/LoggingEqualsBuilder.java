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

package edu.cornell.med.icb.util;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.util.Collection;

/**
 * Like Apache Commons Lang EqualsBuilder but logs if a value becomes non-equals().
 * This will only log the FIRST non-equals() field.
 * This does NOT support the reflection methods of EqualsBuilder.
 *
 * @author Kevin Dorff
 */
public class LoggingEqualsBuilder {

    /**
     * Used to log debug and informational messages.
     */
     private static final Logger LOG = Logger.getLogger(LoggingEqualsBuilder.class);

    /**
     * We use EqualsBuilder under the covers.
     */
    public EqualsBuilder base;

    /**
     * The canonical name of the class we are checking equality for.
     */
    public String classToEquals;

    /**
     * What log level the user wanted non-equals to log to.
     */
    public Level logLevel;

    /**
     * Create a new LoggingEqualsBuilder.
     * @param logLevel The org.apache.log4j.Level to log at non-equals messages at.
     */
    public LoggingEqualsBuilder(
            final Class classToEquals, final Level logLevel) {
        this.logLevel = logLevel;
        base = new EqualsBuilder();
        this.classToEquals = classToEquals.getName();
    }

    public LoggingEqualsBuilder appendSuper(final String fieldName, final boolean superEquals) {
        if (!base.isEquals()) {
            return this;
        }
        base.appendSuper(superEquals);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false");
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final Object a, final Object b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            if (a instanceof Collection) {
                LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                        ".equals() evaluated to false. lhs=" +
                        ArrayUtils.toString(a.toString()) + " rhs=" +
                        ArrayUtils.toString(b.toString()));
            } else {
                LOG.log(logLevel,"Field " + classToEquals + "." + fieldName +
                        ".equals() evaluated to false. lhs=" +
                        a.toString() + " rhs=" + b.toString());
            }
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final long a, final long b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + a + " rhs=" + b);
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final int a, final int b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + a + " rhs=" + b);
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final short a, final short b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + a + " rhs=" + b);
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final char a, final char b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + a + " rhs=" + b);
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final byte a, final byte b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + a + " rhs=" + b);
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final double a, final double b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + a + " rhs=" + b);
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final float a, final float b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + a + " rhs=" + b);
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final boolean a, final boolean b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + a + " rhs=" + b);
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final Object[] a, final Object[] b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + ArrayUtils.toString(a) +
                    " rhs=" + ArrayUtils.toString(a));
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final long[] a, final long[] b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + ArrayUtils.toString(a) +
                    " rhs=" + ArrayUtils.toString(a));
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final int[] a, final int[] b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + ArrayUtils.toString(a) +
                    " rhs=" + ArrayUtils.toString(a));
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final short[] a, final short[] b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + ArrayUtils.toString(a) +
                    " rhs=" + ArrayUtils.toString(a));
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final char[] a, final char[] b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + ArrayUtils.toString(a) +
                    " rhs=" + ArrayUtils.toString(a));
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final byte[] a, final byte[] b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + ArrayUtils.toString(a) +
                    " rhs=" + ArrayUtils.toString(a));
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final double[] a, final double[] b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + ArrayUtils.toString(a) +
                    " rhs=" + ArrayUtils.toString(a));
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final float[] a, final float[] b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + ArrayUtils.toString(a) +
                    " rhs=" + ArrayUtils.toString(a));
        }
        return this;
    }
    
    public LoggingEqualsBuilder append(final String fieldName, final boolean[] a, final boolean[] b) {
        if (!base.isEquals()) {
            return this;
        }
        base.append(a, b);
        if (!base.isEquals()) {
            LOG.log(logLevel, "Field " + classToEquals + "." + fieldName +
                    ".equals() evaluated to false. lhs=" + ArrayUtils.toString(a) +
                    " rhs=" + ArrayUtils.toString(a));
        }
        return this;
    }
    
    public boolean isEquals() {
        return base.isEquals();
    }
}

