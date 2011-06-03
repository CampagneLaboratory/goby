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

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectIterator;

import java.util.List;

/**
 * @author Fabien Campagne
 *         Date: Mar 26, 2011
 *         Time: 7:32:22 PM
 */
public class ColumnFields extends ObjectArraySet<ColumnField> {

    private static final long serialVersionUID = -6258574769357583135L;
    final Object2ObjectOpenHashMap<CharSequence, ColumnField> map =
            new Object2ObjectOpenHashMap<CharSequence, ColumnField>();
    final Int2ObjectOpenHashMap<ColumnField> indexMap =
            new Int2ObjectOpenHashMap<ColumnField>();

    private int previousSize = -1;

    public boolean hasFieldName(CharSequence id) {
        rebuildMap();

        return map.containsKey(id);
    }

    private void rebuildMap() {
        final int setSize = size();
        if (map.size() != setSize) {
            map.clear();
            indexMap.clear();
            for (ColumnField field : this) {
                map.put(field.id, field);
                indexMap.put(field.globalFieldIndex, field);
            }
            previousSize = map.size();
        }
    }

    public ColumnField find(final CharSequence id) {
        rebuildMap();

        return map.get(id);
    }

    public ColumnField find(final int fieldIndex) {
        return indexMap.get(fieldIndex);
    }

    public final ColumnField get(final int index) {
        return list.get(index);
    }

    final ObjectArrayList<ColumnField> list = new ObjectArrayList<ColumnField>();

    @Override
    public ObjectIterator<ColumnField> iterator() {

        rebuildList();
        return list.listIterator();
    }

    public final void rebuildList() {
        final int setSize = size();
        if (previousSize != setSize) {
            list.clear();
            ObjectIterator<ColumnField> it = super.iterator();
            while (it.hasNext()) {
                ColumnField next = it.next();
                list.add(next);
            }
            previousSize = setSize;
        }
    }
}
