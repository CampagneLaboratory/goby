/**
 * Copyright (C) 2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iomanip>
#include <iostream>

#include "goby/Alignments.pb.h"
#include "goby/Alignments.h"

using namespace std;

int main(int argc, const char *const argv[]) {
  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  if (argc != 2) {
    cerr << "usage: " << argv[0] << " <basename>" << endl;
    return -1;
  }

  const string basename = goby::Alignment::getBasename(argv[1]);
  cout << boolalpha << fixed << setprecision(4);

  goby::AlignmentReader alignmentReader = goby::AlignmentReader(basename);

  // write the header lines
  cout << "queryId\treferenceId\treferenceLength\tnumberOfIndels\tnumberOfMismatches\tscore\tstartPosition\talignmentLength\tmatchesReverseStrand" << endl;
  cout << "@HD\tVN:1.0" << endl;
  cout << "@PG\tGoby\tVN:cplusplus" << endl;

  // write the target identifiers and lengths using the target mapping information
  const LIBGOBY_HASH_MAP<string, unsigned>& targetIdentifierToIndex = alignmentReader.getTargetIdentifiers();
  LIBGOBY_HASH_MAP<unsigned, string> targetIndexToIdentifier;
  const vector<unsigned> targetLengths = alignmentReader.getTargetLengths();
  for (LIBGOBY_HASH_MAP<string, unsigned>::const_iterator iter = targetIdentifierToIndex.begin(); iter != targetIdentifierToIndex.end(); ++iter) {
    const string targetId = iter->first;
    const unsigned targetIndex = iter->second;
    targetIndexToIdentifier[targetIndex] = targetId;
    cout << "@SQ\tSN:" << targetId;
    if (!targetLengths.empty()) {
      cout << "\tLN:" << targetLengths.at(targetIndex);
    }
    cout << endl;
  }

  // write the actual alignment information
  for (goby::AlignmentEntryIterator iter = alignmentReader.begin(); iter != alignmentReader.end(); iter++) {
    const goby::AlignmentEntry entry = *iter;
    unsigned targetIndex = entry.target_index();
    int targetLength;
    if (!targetLengths.empty()) {
      targetLength = targetLengths.at(targetIndex);
    } else {
      targetLength = -1;
    }

    cout << entry.query_index()                                      // queryId
      << "\t" << targetIndexToIdentifier.find(targetIndex)->second   // referenceId
      << "\t" << targetLength                                        // referenceLength
      << "\t" << entry.number_of_indels()                            // numberOfIndels
      << "\t" << entry.number_of_mismatches()                        // numberOfMismatches
      << "\t" << entry.score()                                       // score
      << "\t" << entry.position()                                    // startPosition
      << "\t" << entry.query_aligned_length()                        // alignmentLength 
      << "\t" << entry.matching_reverse_strand()                     // matchesReverseStrand
      << endl;
  }

  // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
