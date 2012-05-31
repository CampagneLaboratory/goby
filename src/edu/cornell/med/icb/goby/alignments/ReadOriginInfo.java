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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.objects.Object2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.util.List;

/**
 * Stores and queries ReadOriginInfo instances.
 *
 * @author Fabien Campagne
 *         Date: 4/6/12
 *         Time: 3:13 PM
 */
public class ReadOriginInfo {
    private final List<Alignments.ReadOriginInfo> list;
    private final Object2ObjectAVLTreeMap<Integer, Alignments.ReadOriginInfo> map;

    /**
     * Construct a query object from a list of protocol buffer instances.
     *
     * @param list PB instances (from the alignment header).
     */
    public ReadOriginInfo(final List<Alignments.ReadOriginInfo> list) {
        this.list = list;
        map = new Object2ObjectAVLTreeMap<Integer, Alignments.ReadOriginInfo>();
        for (final Alignments.ReadOriginInfo roi : list) {
            map.put(roi.getOriginIndex(), roi);
        }
    }

    /**
     * Return the protocol buffer ReadOriginInfo list.
     *
     * @return the protocol buffer ReadOriginInfo list.
     */
    public List<Alignments.ReadOriginInfo> getPbList() {
        return list;
    }

    /**
     * Get the read origin info corresponding to the index.
     *
     * @param readOriginIndex index of the object to retrieve.
     * @return the protocol buffer ReadOriginInfo instance corresponding to readOriginIndex.
     */
    public Alignments.ReadOriginInfo getInfo(final int readOriginIndex) {
        return map.get(readOriginIndex);
    }

    /**
     * Return a list of protocol buffer ReadOriginInfo builders.
     * @return list of protocol buffer ReadOriginInfo builders.
     */
    public ObjectArrayList<Alignments.ReadOriginInfo.Builder> getPBBuilderList() {
        ObjectArrayList<Alignments.ReadOriginInfo.Builder> result = new ObjectArrayList<Alignments.ReadOriginInfo.Builder>();
        for (Alignments.ReadOriginInfo e : list) {
            result.add(Alignments.ReadOriginInfo.newBuilder(e));
        }
        return result;
    }

    /**
     * The number of read origin info/read groups defined.
     *
     * @return 0 if no read origin info are defined.
     */
    public int size() {
        return list.size();
    }
}
