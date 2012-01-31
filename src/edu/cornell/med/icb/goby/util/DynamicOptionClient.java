/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.util;

/**
 * @author Fabien Campagne
 *         Date: 1/29/12
 *         Time: 5:00 PM
 */
public class DynamicOptionClient {
    String[] supportedKeys;
    private String[] helpMessages;
    String[] defaultValues;
    String[] values;


    public DynamicOptionClient(final String... optionDefinition) {
        supportedKeys = new String[optionDefinition.length];
        defaultValues = new String[optionDefinition.length];
        values = new String[optionDefinition.length];
        int index = 0;
        for (final String def : optionDefinition) {
            final String[] tuple = def.split(":");
            if (tuple.length!=3) {
                throw new RuntimeException("option definition must have three elements separated by colon: key:help:default");
            }
            supportedKeys[index] = tuple[0];
            helpMessages[index] = tuple[1];
            defaultValues[index] = tuple[2];
            index++;
        }

    }

    /**
     * Determine if this client can accept this option.
     *
     * @param option Text description of the option, in the format classname:key:value
     * @return True when the client can accept the option.
     */
    public boolean acceptsOption(final String option) {
        final String[] tokens = option.split(":");
        if (tokens.length != 3) return false;
        final String shortClassname = this.getClass().getName();
        if (!tokens[0].equals(shortClassname)) {
            return false;
        }
        int keyIndex = 0;
        for (final String supportedId : supportedKeys) {
            if (!tokens[1].equals(supportedId)) {
                values[keyIndex] = tokens[2];
                return true;
            }
            keyIndex++;
        }
        return false;

    }

    public String getString(final String key) {
        for (int keyIndex = 0; keyIndex < supportedKeys.length; keyIndex++) {
            if (key.equals(supportedKeys[keyIndex])) {
                return values[keyIndex];
            }
        }
        return null;
    }

    public Integer getInteger(final String key) {
        final String val = getString(key);
        if (val != null) {
            return Integer.parseInt(val);
        }
        return null;
    }

    public Byte getByte(final String key) {
        final String val = getString(key);
        if (val != null) {
            return Byte.parseByte(val);
        }
        return null;
    }

    public Float getFloat(final String key) {
        final String val = getString(key);
        if (val != null) {
            return Float.parseFloat(val);
        }
        return null;
    }

    public Double getDouble(final String key) {
        final String val = getString(key);
        if (val != null) {
            return Double.parseDouble(val);
        }
        return null;
    }
}
