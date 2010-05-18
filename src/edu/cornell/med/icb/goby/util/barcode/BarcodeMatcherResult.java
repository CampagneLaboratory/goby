/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.util.barcode;

import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;

/**
 * A match from a BarcodeMatcher.
 *
 * @author Kevin Dorff
 */
public class BarcodeMatcherResult {
    /** Logging. */
    private static final Log LOG = LogFactory.getLog(BarcodeMatcherResult.class);

    private int barcodeIndex = -1;
    private int numberOfDiffs = -1;
    private int sequenceStartPosition = -1;
    private int sequenceLength = -1;
    private int barcodeStartPosition = -1;
    private int barcodeMatchLength = -1;
    private boolean ambiguous = false;

    public BarcodeMatcherResult(
            final int barcodeIndex,
            final int numberOfDiffs,
            final int sequenceStartPosition, final int sequenceLength,
            final int barcodeStartPosition, final int barcodeMatchLength) {
        this.barcodeIndex = barcodeIndex;
        this.numberOfDiffs = numberOfDiffs;
        this.sequenceStartPosition = sequenceStartPosition;
        this.sequenceLength = sequenceLength;
        this.barcodeStartPosition = barcodeStartPosition;
        this.barcodeMatchLength = barcodeMatchLength;
    }

    private static Map<String, Field> FIELDS_MAP = null;

    public static synchronized void ensureFieldsMap() {
        if (FIELDS_MAP != null) {
            return;
        }
        final Field[] fields = BarcodeMatcherResult.class.getDeclaredFields();
        FIELDS_MAP = new HashMap<String, Field>(fields.length);
        for (final Field field : fields) {
            FIELDS_MAP.put(field.getName(), field);
        }
    }

    /**
     * This constructor is provided as an alternate to the one above. It uses a map
     * to set the values of the fields. This method is considerably slower because it uses reflection to
     * set the values of the class from a map so don't use this in the case where you are creating
     * millions of these objects.
     * @param init the map to initiaize the class with
     */
    public BarcodeMatcherResult(final Map<String, Object> init) {
        ensureFieldsMap();
        for (final Map.Entry<String, Object> entry : init.entrySet()) {
            try {
                final Object value = entry.getValue();
                if (value == null) {
                    continue;
                }
                final Field f = FIELDS_MAP.get(entry.getKey());
                if (value instanceof Integer) {
                    f.setInt(this, (Integer) value);
                } else if (value instanceof Boolean) {
                    f.setBoolean(this, (Boolean) value);
                }
            } catch (IllegalAccessException e) {
                LOG.error(e);
            } catch (IllegalArgumentException  e) {
                LOG.error(e);
            }
        }
    }

    public int getBarcodeIndex() {
        return barcodeIndex;
    }

    public int getNumberOfDiffs() {
        return numberOfDiffs;
    }

    public int getSequenceStartPosition() {
        return sequenceStartPosition;
    }

    public int getSequenceLength() {
        return sequenceLength;
    }

    public int getBarcodeStartPosition() {
        return barcodeStartPosition;
    }

    public int getBarcodeMatchLength() {
        return barcodeMatchLength;
    }

    public boolean isAmbiguous() {
        return ambiguous;
    }

    public void setAmbiguous(final boolean ambiguous) {
        this.ambiguous = ambiguous;
    }

    public MutableString sequenceOf(final MutableString sequence) {
        return sequence.substring(
            getSequenceStartPosition(), getSequenceStartPosition() + getSequenceLength());
    }

    public MutableString sequenceBarcodeOf(final MutableString sequence) {
        return sequence.substring(
            getBarcodeStartPosition(), getBarcodeStartPosition() + getBarcodeMatchLength());
    }

    public MutableString barcodeWithMatchingAdapterOf(final BarcodeMatcher matcher) {
        return matcher.getBarcodeWithAdapterAtIndex(getBarcodeIndex()).substring(0, getBarcodeMatchLength());
    }

    public MutableString barcodeOnlyOf(final BarcodeMatcher matcher) {
        return matcher.getBarcodeOnlyAtIndex(getBarcodeIndex());
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final BarcodeMatcherResult that = (BarcodeMatcherResult) o;

        if (ambiguous != that.ambiguous) {
            return false;
        }
        if (barcodeIndex != that.barcodeIndex) {
            return false;
        }
        if (barcodeMatchLength != that.barcodeMatchLength) {
            return false;
        }
        if (barcodeStartPosition != that.barcodeStartPosition) {
            return false;
        }
        if (numberOfDiffs != that.numberOfDiffs) {
            return false;
        }
        if (sequenceLength != that.sequenceLength) {
            return false;
        }
        if (sequenceStartPosition != that.sequenceStartPosition) {
            return false;
        }

        return true;
    }

    @Override
    public int hashCode() {
        int result = barcodeIndex;
        result = 31 * result + numberOfDiffs;
        result = 31 * result + sequenceStartPosition;
        result = 31 * result + sequenceLength;
        result = 31 * result + barcodeStartPosition;
        result = 31 * result + barcodeMatchLength;
        result = 31 * result + (ambiguous ? 1 : 0);
        return result;
    }
}
