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

package edu.cornell.med.icb.goby.readers.vcf;

import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;

/**
 * Stores association from columns to groups.
 *
 * @author Fabien Campagne
 *         Date: 5/20/12
 *         Time: 12:32 PM
 */
public class GroupAssociations {
    private ObjectArrayList<String> associations;

    private Object2ObjectMap<String, ObjectArraySet<String>> groupToColumns = new Object2ObjectOpenHashMap<String, ObjectArraySet<String>>();

    public GroupAssociations(String associationsAsText, ColumnInfo formatColumn, String[] sampleIds) {
       // parse associations stored in text in the header:
        associations = parse(associationsAsText);
        if (formatColumn != null) {
        // now associate each sampleId to the columns that result from the combination of FORMAT fields and sample ids:
        for (String sample : sampleIds) {
            for (final ColumnField formatField : formatColumn.fields) {
                associate(String.format("%s[%s]", sample, formatField.id), sample);
            }
        }
        }
    }

    private ObjectArrayList<String> parse(String associationsAsText) {
        ObjectArrayList<String> associations = new ObjectArrayList<String>();
        if (associationsAsText == null) return associations;
        String[] tokens = associationsAsText.split(",");
        for (int i = 0; i < tokens.length; i++) {
            String token = tokens[i];
            associations.add(token);
            String[] col_group = token.split("=");
            associate(col_group[0], col_group[1]);
        }
        return associations;
    }

    private void associate(String columnName, String group) {
        ObjectArraySet<String> colList = groupToColumns.get(group);
        if (colList == null) {
            colList = new ObjectArraySet<String>();
            groupToColumns.put(group, colList);
        }
        colList.add(columnName);
    }

    public ObjectArraySet<String> getColumnsWithGroup(String group) {
        return groupToColumns.get(group);
    }

    /**
     * Return a comma separated list of groups associated with a column.
     *
     * @param columnName Name of the column.
     * @return comma separated list of group names.
     */
    public String listGroupsAsString(final String columnName) {
        final StringBuffer sb = new StringBuffer();
        boolean first = true;
        for (final String association : associations) {

            if (isMatchingColumnName(columnName, association)) {
                if (!first) {
                    sb.append(",");
                }
                sb.append(association.replace(getMatchingString(columnName, association), ""));
                first = false;
            }
        }
        return sb.toString();
    }

    private String getMatchingString(String columnName, String association) {
        if (association.startsWith(columnName)) {
            return columnName + '=';
        } else if (association.startsWith("INFO/") && association.contains(columnName)) {
            return "INFO/" + columnName + '=';
        }
        String strippedColName = columnName.substring("INFO[".length(),columnName.length()-1);
        if (association.startsWith("INFO") && association.contains(strippedColName)) {
            return "INFO/" + strippedColName + '=';// association.substring(0,association.indexOf('=')+1) ;
        } else return "not-matchjing=q3p-ow";

    }

    private boolean isMatchingColumnName(String columnName, String association) {
        if (association==null || columnName==null) {
            return false;
        }
        final String infoColumn = "INFO/" + columnName + "=";

        if (association.startsWith(columnName) || association.startsWith(infoColumn)) return true;
        String strippedColName = columnName.replace("INFO[", "").replace("]", "");

        return association.startsWith("INFO") && association.contains(strippedColName);
    }

    /**
     * Return a list of groups associated with a column.
     *
     * @param columnName Name of the column.
     * @return list of group names that column is associated with.
     */
    public ObjectArrayList<String> listGroups(final String columnName) {
        final ObjectArrayList<String> result = new ObjectArrayList<String>();
        for (final String association : associations) {

            if (isMatchingColumnName(columnName, association)) {

                final String matchingString = getMatchingString(columnName, association);
                result.add(association.replace(matchingString, ""));
            }

        }
        return result;
    }


}
