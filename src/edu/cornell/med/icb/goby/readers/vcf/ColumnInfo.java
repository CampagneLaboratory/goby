/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.readers.vcf;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;

import java.util.Arrays;

/**
 * Meta-information about columns of a VCF file. This information is parsed from the VCF header.
 *
 * @author Fabien Campagne
 *         Date: Mar 26, 2011
 *         Time: 7:31:33 PM
 */
public class ColumnInfo implements Cloneable {
    public String getColumnName() {
        return columnName;
    }

    /**
     * Name of the column.
     */
    protected String columnName;
    /**
     * Set of fields for the column.
     */
    public ColumnFields fields = new ColumnFields();
    /**
     * Index of this column in the file being parsed. Index is zero-based. Zero is the first (left-most) column.
     * Value is -1 if this index has not been set yet.
     */
    public int columnIndex = -1;
    public boolean useFormat = false;
    public int formatIndex = 0;

    /**
     * Create a ColumnInfo with provided information. This method does not transfer groups from field to column.
     *
     * @param columnName Name of the new column.
     * @param fields     Fields in this column.
     */
    public ColumnInfo(String columnName, ColumnField... fields) {
    this(columnName,false,fields);
    }
    /**
     * Create a ColumnInfo with provided information.
     *
     * @param columnName Name of the new column.
     * @param transferGroups  Indicate that the column groups should be initialized with a copy of the fields' groups.
     * @param fields     Fields in this column.
     */
    public ColumnInfo(String columnName, final boolean transferGroups, ColumnField... fields) {
        this.columnName = columnName;
        this.fields.addAll(Arrays.asList(fields));
        for (ColumnField field : fields) {
            field.column = this;
            if (transferGroups) {
                columnGroups.addAll(field.getGroups());
            }
        }

    }

    private ObjectArraySet<String> columnGroups = new ObjectArraySet<String>();


    /**
     * Create a ColumnInfo with empty information.
     */
    public ColumnInfo() {

    }

    /**
     * Add a field to this column. Note that any group associated with the column is added to the field.
     * @param field Field to add to this column.
     */
    public void addField(ColumnField field) {
        fields.add(field);
        field.column = this;
        for (String columnGroup: columnGroups) {
            field.addGroup(columnGroup);
        }
    }

    public boolean hasField(String fieldName) {
        return fields.hasFieldName(fieldName);
    }

    public ColumnField getField(String fieldName) {
        return fields.find(fieldName);
    }


    public ColumnInfo copy() {
        ColumnField[] fieldCopy = new ColumnField[fields.size()];
        int i = 0;
        for (ColumnField f : fields) {
            fieldCopy[i++] = new ColumnField(f.id, f.numberOfValues, f.type, f.description,
                    f.getGroups().toArray(new String[f.getGroups().size()]));
        }
        ColumnInfo copy = new ColumnInfo(columnName, fieldCopy);
        copy.addGroup(columnGroups.toArray(new String[columnGroups.size()]));
        return copy;
    }

    private void addGroup(String... columnGroups) {
       this. columnGroups.addAll(ObjectArrayList.wrap(columnGroups));
    }

    public String[] getGroups() {
        return columnGroups.toArray(new String[columnGroups.size()]);
    }
}
