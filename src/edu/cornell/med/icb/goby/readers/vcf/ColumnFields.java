/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.readers.vcf;

import it.unimi.dsi.fastutil.objects.ObjectArraySet;

/**
 * @author Fabien Campagne
*         Date: Mar 26, 2011
*         Time: 7:32:22 PM
*/
public class ColumnFields extends ObjectArraySet<ColumnField> {
    public boolean hasFieldName(CharSequence id) {
        for (ColumnField field : this) {
            if (id.equals(field.id))
                return true;
        }
        return false;
    }
    public ColumnField find(CharSequence id) {
        for (ColumnField field : this) {
            if (id.equals(field.id))
                return field;
        }
        return null;
    }
}
