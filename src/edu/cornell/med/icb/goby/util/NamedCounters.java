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

import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.*;

/**
 * @author Fabien Campagne
 *         Date: 1/13/12
 *         Time: 2:36 PM
 */
public class NamedCounters {
    private ObjectArrayList<LongNamedCounter[]> counters;
    int numFiles;

    private ObjectSet<String> names = new ObjectOpenHashSet<String>();
    private Object2IntMap<String> nameToIndex = new Object2IntOpenHashMap<String>();

    @Override
    public String toString() {
        StringBuffer buffer = new StringBuffer();

        for (String name : names) {
            buffer.append(String.format("[%s: %s]%n", name, ObjectArrayList.wrap( getArray(name))));
        }
        return buffer.toString();
    }

    public NamedCounters(int numFiles) {

        this.counters = new ObjectArrayList<LongNamedCounter[]>();
        this.numFiles = numFiles;
    }

    private int index = 0;

    public void register(String name, int numFiles) {
        if (!names.contains(name)) {
            counters.add(new LongNamedCounter[numFiles]);

            nameToIndex.put(name, index);
            for (int fileIndex = 0; fileIndex < numFiles; fileIndex++) {
                counters.get(index)[fileIndex] = new LongNamedCounter(name);
            }
            index++;
        }
    }

    /**
     * Returns an array of LongNamedCounter, indexed by fileIndex.
     *
     * @param name
     * @return
     */
    public LongNamedCounter[] getArray(String name) {
        return counters.get(nameToIndex.getInt(name));
    }

    public LongNamedCounter get(String name, int fileIndex) {
        assert nameToIndex.containsKey(name) :"counter name must exist: "+name ;
        int counterIndex = nameToIndex.getInt(name);
        return counters.get(counterIndex)[fileIndex];
    }

    public String[] ids() {
        return nameToIndex.keySet().toArray(new String[nameToIndex.size()]);
    }
}
