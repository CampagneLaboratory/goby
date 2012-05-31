#!/bin/env groovy

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

//
// When files are sequenced, the directory containing the sequences (s_#_sequence.txt.gz)
// contains a file "Summary.xml". This will parse that file and display the
// sampleId's and sequenceId's within the summary file.
//
// As command line arguments, this script takes one or more summary files
// or directories. Any directories that are provided will be used to find
// the file Summary.xml.
//

import org.xml.sax.SAXParseException

for (arg in args) {
    def file = new File(arg)
    if (!file.exists()) {
        println "Argument ${arg} doesn't represent either a file or a directory"
        continue
    }
    if (file.isDirectory()) {
        file = new File(file, "Summary.xml")
    }
    parseXmlSummaryFile(file)
}


def parseXmlSummaryFile(File summaryFile) {
    if (!summaryFile.exists()) {
        println "Could not find file to parse: ${summaryFile.toString()}"
        return
    }
    println "Summary file: ${summaryFile.toString()}"
    try {
        def xml = new XmlParser().parse(summaryFile)
        def lanes = xml.Samples.Lane
        for (lane in lanes) {
            String laneNumber = lane.laneNumber.text().trim()
            String sampleId = lane.sampleId.text().trim()
            println "   sampldId = ${sampleId},   laneNumber=${laneNumber}"
        }
    } catch (SAXParseException e) {
        println "WARN: Couldn't parse ${summaryFile.toString()}"
    }

}