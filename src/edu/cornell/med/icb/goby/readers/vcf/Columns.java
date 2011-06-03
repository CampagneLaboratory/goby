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

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectIterator;
import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: Mar 26, 2011
 *         Time: 7:33:07 PM
 */
public class Columns extends ObjectArraySet<ColumnInfo> {
    private int previousSize;

    public boolean hasColumnName(CharSequence id) {
        for (ColumnInfo info : this) {
            if (id.equals(info.columnName))
                return true;
        }
        return false;
    }

    public ColumnInfo find(String columnName) {
        for (ColumnInfo info : this) {
            if (columnName.equals(info.columnName))
                return info;
        }
        return null;
    }

    final ObjectArrayList<ColumnInfo> list = new ObjectArrayList<ColumnInfo>();

    @Override
    public final ObjectIterator<ColumnInfo> iterator() {
        rebuildList();
        return list.listIterator();
    }

    private void rebuildList() {
        final int setSize = size();
        if (previousSize != setSize) {
            list.clear();
            ObjectIterator<ColumnInfo> it = super.iterator();
            while (it.hasNext()) {
                ColumnInfo next = it.next();
                list.add(next);
            }
            previousSize = setSize;
        }
    }
}
