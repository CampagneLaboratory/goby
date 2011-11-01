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

/**
 * Meta information about fields of VCF columns.
 *
 * @author Fabien Campagne
 *         Date: Mar 26, 2011
 *         Time: 3:17:32 PM
 */
public class ColumnField {
    /**
     * Field identifier.
     */
    public String id;
    /**
     * Global index of this field across all columns of the input file.
     * Value -1 indicates that the index has not yet been assigned.
     */
    public int globalFieldIndex=-1;
    /**
     * Number of values listed in this field.
     */
    public int numberOfValues;
    /**
     * Data type for the column.
     */
    public ColumnType type;
    /**
     * Description for this field of the column.
     */
    public String description;
    /**
     * An optional group to which this column field belongs. The group is supported so that user interfaces can
     * manipulate columns that are related together. For instance, a UI could provide the option of hidding or
     * showing all columns in a group at the same time. The default group is "MAIN".
     */
    public String group="MAIN";

    
    public ColumnInfo column;
    public boolean notSeen;

    /**
     * Construct an empty field.
     */
    public ColumnField() {
    }

    /**
     * Construct a field with provided information.
     * @param id Field identifier
     * @param numberOfValues Number of values in the field
     * @param type  Field type.
     * @param description Description for field content.
     */
    public ColumnField(String id, int numberOfValues, ColumnType type, String description) {
        this.id = id;
        this.numberOfValues = numberOfValues;
        this.type = type;
        this.description = description;
    }

}
