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
// After files have been sequenced by one or more Flow Cells
// (and de-multiplexed) and transferred from the facility that
// does the sequencing and converted to .compact-reads files,
// you are left with a directory structure similar to
//    + FlowCellA
//      + 001
//        + s_1_sequence.compact-reads
//        + ..
//        + s_8_sequence.compact-reads
//        + Summary.xml
//      + 002
//        + s_1_sequence.compact-reads
//        + ..
//        + s_8_sequence.compact-reads
//        + Summary.xml
//    + FlowCellB
//      + 001
//        + s_1_sequence.compact-reads
//        + ..
//        + s_8_sequence.compact-reads
//        + Summary.xml
//      + 002
//        + s_1_sequence.compact-reads
//        + ..
//        + s_8_sequence.compact-reads
//        + Summary.xml
//
// Combine this file structure with a tsv file which describes the
// data that has been sequenced, the --project-tsv file. This file
// The Summary.xml file converts from lane id (such as s_3_sequence.compact-reads
// comes from lane 3) to the sequencer-id (such as 32, which is generally something
// like S32 in the actual project). The --project-tsv file should have a column
// which is the sequencer-id, you specify which column this is with the
// --sample-id-column argument. Specify the columns you want to rename the file
// to with the --rename-columns argument. An example command for this script is,
// running in dry-run mode (no renames will take place):
//
//     mkdir output-dir
//     ./rename.groovy --dry-run -o ./output_dir -i sequencerId -p project-data.tsv
//         -r shortSampleIdentifier typeOfSample phenotype _lane_ %lane _fc_ %flowCell
//         -d FlowCell?/0*
//
// Run with --help to get complete usage information.
//

import org.apache.commons.cli.Option

class RenameSequencedFiles {

    def FLOW_CELL_PATTERN = ~/Flow_*Cell_*([a-zA-Z]+)/
    // The following PROCESS_SUFFIX should NOT be a regex
    def PROCESS_SUFFIX = "compact-reads" // "txt.gz"
    def READS_FILE_PATTERN = ~/^s_(\d+)_sequence\.${PROCESS_SUFFIX}$/

    boolean dryRun
    File outputFolderFile
    File projectTsvFile
    String sampleIdColumn
    List<String> renameColumns
    List<ProcessDetails> processDetailsList
    Map<String, List<String>> sequencerIdToRenameColumns = [:]
    // Special columns names are data that come from filenames or directory names
    List<String> specialColumnNames = ["%flowCell", "%lane"]

    public static void main(final String[] args) {
        new RenameSequencedFiles().execute(args)
    }

    void execute(String[] args) {
        // Configure command line args, parse --project-tsv, read Summary.xml files
        if (!configure(args)) {
            return
        }
        // Given the configuration, figure out the rename from's and to's
        Map<String, String> renameFromToMap = createRenameMap()
        if (!renameFromToMap) {
            return
        }
        if (this.dryRun) {
            // Dry run output
            renameFromToMap.each { String k, String v ->
                println "'${k}'  would be renamed to  '${v}'"
            }
        } else {
            def renameEntries = renameFromToMap.entrySet()
            // Verify the files to rename, source exists and dest doesn't exist
            for (renameEntry in renameEntries) {
                File source = new File(renameEntry.key)
                File dest = new File(renameEntry.value)
                if (!source.exists()) {
                    System.err.println "The file to rename FROM for ${source.toString()} to ${dest.toString()}  doesn't exist."
                    return
                } else if (dest.exists()) {
                    System.err.println "The file to rename TO for ${source.toString()} to ${dest.toString()}  already exist."
                    return
                }
            }
            // Perform the rename, then verify that it worked, the source file no longer exists and the
            // dest file now does exist.
            for (renameEntry in renameEntries) {
                println "Moving/Renaming:  '${renameEntry.key}'  to  '${renameEntry.value}'"
                File source = new File(renameEntry.key)
                File dest = new File(renameEntry.value)
                def renamed = source.renameTo(dest)
                if (!renamed) {
                    System.err.println "Move/Rename failed. It is likely you tried to move/rename across different filesystems."
                    return
                }
                if (source.exists()) {
                    System.err.println "The file to rename FROM exists after the rename, which it shouldn't."
                    return
                } else if (!dest.exists()) {
                    System.err.println "The file to rename TO doesn't exist after the rename, which it should."
                    return
                }
            }
        }
    }
    
    boolean configure(final String[] args) {
        // Command line parsing configuration
        def cli = new CliBuilder(usage: 'rename.groovy')
        cli.h(longOpt: "help", "help")
        cli.y(longOpt: "dry-run", "Execute 'dry-run' (don't actually rename / move files).")
        cli.o(longOpt: "output-folder", args:1, required: true, "The output folder, which must be on the filesystem as the source files.")
        cli.p(longOpt: "project-tsv", args:1, required: true, "The input tsv file that contains the project data.")
        cli.i(longOpt: "sample-id-column", args:1, required: true,
            "The name of the column that should exactly match the value " +
            "from Summery/Samples/Lane/sampleId from the Summary.xml file.")
        cli.r(longOpt: "rename-columns", 
            args: Option.UNLIMITED_VALUES, required: true, valueSeparator: ' ' as char,
            "The names of the columns from the header line of the --project-tsv file, " +
            "separated by spaces, used to construct the new filenames. " +
            "String literals can be inlcuded by flanking the string with '_' and using '_' for " +
            "spaces, such as '_cat_and_dog_' would be the literal 'cat and dog'. " +
            "Special column names are ${specialColumnNames.join(', ')}. " +
            "If you specify a single value @columns.txt, the columns can be read " +
            "file the named file, one column per line.")
        cli.d(longOpt: "process-dirs", 
            args: Option.UNLIMITED_VALUES, required: true, valueSeparator: ' ' as char,
            "The directories to process, each specified directory should contain " +
            "one or more .${PROCESS_SUFFIX} files, named in the " +
            "format of s_#_sequence.${PROCESS_SUFFIX}, " +
            "where # specifies the lane number, and a Summary.xml file. " +
            "Each specified directory should contain the Flow Cell name, such as 'FlowCellA/006' so " +
            "the script can determine from which Flow Cell the reads in each directory come. " +
            "If you specify a single value @directories.txt, the directories can be read " +
            "file the named file, one directory per line, but this method does not support wildcards.")
        def options = cli.parse(args)
        if (!options) {
            // Options didn't parse
            return false
        }
        if (options.h) {
            cli.usage()
            return false
        }
        // Retrieve pared command line options
        this.dryRun = options.y
        this.outputFolderFile = new File(options.o)
        if (!this.outputFolderFile.exists() || !this.outputFolderFile.isDirectory()) {
            System.err.println "Directory --output-folder doesn't exist or it isn't a directory or it isn't writable."
            return false
        }
        this.projectTsvFile = new File(options.p)
        this.sampleIdColumn = options.i
        this.renameColumns = options.rs
        this.processDetailsList = []
        def directories = options.ds
        // Verify the directories that were specified to be processed...
        // that they are directories, they contain files, and each directory
        // contains a Summary.xml file.
        for (String directory in directories) {
            File directoryFile = new File(directory)
            if (!directoryFile.exists() || !directoryFile.isDirectory()) {
                System.err.println "Cannot add directory ${directory} because it doesn't exist or isn't a directory."
                return false
            }
            def flowCellName = flowCellFromDirName(directory)
            if (!flowCellName) {
                System.err.println "Could not find the flow cell specified for directory ${directoryFile.toString()}."
                return false
            }
            File summaryFile = new File(directoryFile, "Summary.xml")
            if (!summaryFile.exists() || !summaryFile.isFile()) {
                System.err.println "Cannot find Summary.xml within ${directoryFile.toString()} because it doesn't exist or isn't a file."
                return false
            }
            // Make sure the directory in question has some files that have the appropriate filenames
            List<File> files = filesForDir(directoryFile)
            if (!files) {
                System.err.println "Cannot add directory ${directoryFile.toString()} because doesn't seem to contain any s_#_sequence.${PROCESS_SUFFIX} files."
                return false
            }
            // Parse the XML summary file for a directory of files
            Map<String, String> laneToSequencerIdMap = parseXmlSummaryFile(summaryFile)
            if (!laneToSequencerIdMap) {
                System.err.println "Error parsing Summary.xml file  ${summaryFile.toString()}."
                return false
            }
            processDetailsList << new ProcessDetails('flowCell': flowCellName, 'summaryFile': summaryFile, 'laneToSequencerIdMap': laneToSequencerIdMap, 'files': files)
        }
        if (!processDetailsList) {
            System.err.println "No directories to process."
            return false
        }
        // Make sure the output folder exists and is a directory
        if (!this.outputFolderFile.exists() || !this.outputFolderFile.isDirectory()) {
            System.err.println "The --output-folder ${this.outputFolderFile.toString()} directory does not exist."
            return false
        }
        
        // Some sanity checks on the --project-tsv file exists and is a file
        if (!this.projectTsvFile.exists() || !this.projectTsvFile.isFile()) {
            System.err.println "The --project-tsv ${this.projectTsvFile.toString()} file did not exist or is not a file."
            return false
        }
        if (!this.projectTsvFile.toString().endsWith(".tsv")) {
            System.err.println "The --project-tsv ${this.projectTsvFile.toString()} filename didn't end with .tsv"
            return false
        }
        // Read the project tsv file
        if (!readProjectData()) {
            return false
        }
        return true
    }

    boolean readProjectData() {
        Map<String, Integer> columnNameToNumberMap = [:]
        List<String> projectColumns
        int lineNum = 1
        Integer sampleIdColumnNum
        for (String line in projectTsvFile.readLines()) {
            if (line.startsWith("#")) {
                continue
            }
            if (lineNum == 1) {
                projectColumns = line.split("\t") as List
                // Make a map of column name to column number
                int colNo = 0
                for (String projectColumn in projectColumns) {
                    if (!projectColumn) {
                        System.err.println "The --project-tsv file ${projectTsvFile.toString()} contains a column header that is blank."
                        return false
                    }
                    columnNameToNumberMap[projectColumn] = colNo++
                }
                // Find the sampleIdColumn
                sampleIdColumnNum = columnNameToNumberMap[sampleIdColumn]
                if (sampleIdColumnNum == null) {
                    System.err.println "The --sample-id-column ${sampleIdColumn} not found in --project-tsv ${this.projectTsvFile.toString()} file."
                    return false
                }
                // Verify that the rename columns are contained within projectColumns
                for (String renameColumn in this.renameColumns) {
                    if (renameColumn.startsWith("_") && renameColumn.endsWith("_")) {
                        // Plain text, it's fine.
                    } else if (this.specialColumnNames.contains(renameColumn)) {
                        // Special column name it's fine
                    } else {
                        if (!projectColumns.contains(renameColumn)) {
                            System.err.println "The --rename-columns column ${renameColumn} not found in --project-tsv ${this.projectTsvFile.toString()} file."
                            return false
                        }
                    }
                }
            } else {
                // Verify matching number of columns in data line
                String[] lineParts = line.split("\t")
                if (lineParts.size() != projectColumns.size()) {
                    System.err.println "The --project-tsv file ${projectTsvFile.toString()} line ${lineNum} should contain ${projectColumns.size()} columns but it contains ${lineParts.size()}."
                    return false
                }
                // Determine the rename values for one file
                List<String> renameValues = []
                for (String renameColumn in this.renameColumns) {
                    if (renameColumn.startsWith("_") && renameColumn.endsWith("_")) {
                        renameValues << renameColumn.replaceAll("_", " ").trim()
                    } else if (this.specialColumnNames.contains(renameColumn)) {
                        renameValues << renameColumn
                    } else {
                        renameValues << lineParts[columnNameToNumberMap[renameColumn]].trim()
                    }
                }
                sequencerIdToRenameColumns[lineParts[sampleIdColumnNum]] = renameValues
            }
            lineNum++
        }
        return true
    }

    /**
     * Takes a String or a File and looks for the Flow Cell pattern, such that
     *    "FlowCellA/other"  => "A"
     *    "stuffFlow_Cell__b/other"  => "B"
     *    "stuff/Flow_CellBG/other"  => "BG"
     */
    String flowCellFromDirName(dir) {
        def flowCellName = null
        dir.toString().find(FLOW_CELL_PATTERN) { fullMatch, foundFlowCell ->
            flowCellName = foundFlowCell.toUpperCase()
        }    
        return flowCellName
    }

    /**
     * Given a directory, find the files that match READS_FILE_PATTERN
     */
    List<File> filesForDir(File dirFile) {
        List<File> files = []
        if (dirFile.exists()) {
            dirFile.eachFileMatch(READS_FILE_PATTERN) { File match ->
                files << match
            }
        }
        return files
    }

    /**
     * Given a file in the format s_6_sequence.compact-reads parse out the
     * lane number (6, in this case).
     */
    String laneFromFilename(File file) {
        String filename = file.toString().replaceAll("\\\\","/").split("/")[-1]
        String lane = null
        filename.find(READS_FILE_PATTERN) { String fullmatch, String foundLane ->
            lane = foundLane
        }        
        return lane
    }

    /**
     * Parse a map of sequencer-id to lane number from the
     * Summary.xml file.
     */
    Map<String, String> parseXmlSummaryFile(File summaryFile) {
        Map<String, String> laneToSequencerIdMap = [:]
        def xml = new XmlParser().parse(summaryFile)
        def lanes = xml.Samples.Lane
        for (lane in lanes) {
            String laneNumber = lane.laneNumber.text().trim()
            String sequencerId = lane.sampleId.text().trim()
            laneToSequencerIdMap[laneNumber] = sequencerId
        }            
        return laneToSequencerIdMap
    }

    /**
     * Given the existing configuration and directories,
     * create a rename from-to map.
     */
    Map<String, String> createRenameMap() {
        Map<String, String> renameFromToMap = [:]
        for (ProcessDetails details in processDetailsList) {
            String flowCell = details.flowCell
            Map<String, String> laneToSequencerIdMap = details.laneToSequencerIdMap
            File summaryFile = details.summaryFile
            List<File> files = details.files
            for (file in files) {
                String lane = laneFromFilename(file)
                String sequencerId = laneToSequencerIdMap[lane]
                if (!sequencerId) {
                    System.err.println "Couldn't determine sequencer-id for lane ${lane} for file ${file} from Summary.xml file ${summaryFile}"
                    return null
                }
                List<String> renameValues = []
                def renameColumns = sequencerIdToRenameColumns[sequencerId]
                if (renameColumns != null) {
                    for (String renameColumn in sequencerIdToRenameColumns[sequencerId]) {
                        switch (renameColumn) {
                            case "%flowCell":
                                renameValues << flowCell
                                break
                            case "%lane":
                                renameValues << lane
                                break
                            default:
                                renameValues << renameColumn
                                break
                        }
                    }
                    String oldFilename = file.toString()
                    String newFilename = "${this.outputFolderFile}/${renameValues.join('-')}.${PROCESS_SUFFIX}"
                    if (renameFromToMap.containsKey(oldFilename) ||
                        renameFromToMap.containsValue(oldFilename) ||
                        renameFromToMap.containsKey(newFilename) ||
                        renameFromToMap.containsValue(newFilename)) {
                        System.err.println "Tried to add rename ${oldFilename} to ${newFilename} but one of the filenames was not unique."
                        return null
                    }
                    renameFromToMap[oldFilename] = newFilename
                } else {
                    System.err.println "WARNING: sequencer-id ${sequencerId} not found in --project-tsv file ${projectTsvFile.toString()} so file ${file.toString()} cannot be renamed"
                }
            } 
        }
        return renameFromToMap
    }
}

class ProcessDetails {
    String flowCell
    File summaryFile
    Map<String, String> laneToSequencerIdMap
    List<File> files
}
