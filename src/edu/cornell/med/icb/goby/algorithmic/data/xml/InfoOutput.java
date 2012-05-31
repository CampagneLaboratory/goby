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

package edu.cornell.med.icb.goby.algorithmic.data.xml;


import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import java.util.ArrayList;

/**
 * @author Fabien Campagne
 *         Date: 10/23/11
 *         Time: 5:43 PM
 */
@XmlRootElement(name = "info-output")
public class InfoOutput {

    @XmlElement(name = "al")
    @XmlElementWrapper(name = "annotation-lengths")
    public ArrayList<AnnotationLength> lengths = new ArrayList<AnnotationLength>();

    @XmlElement(name = "tc")
    @XmlElementWrapper(name = "total-counts")
    public ArrayList<SampleTotalCount> totalCounts = new ArrayList<SampleTotalCount>();
}
