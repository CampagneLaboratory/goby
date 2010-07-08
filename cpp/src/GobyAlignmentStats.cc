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

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <locale>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

#include <sys/types.h>
#include <sys/stat.h>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#include "goby/common.h"
#include "goby/Alignments.pb.h"
#include "goby/Alignments.h"
#include "goby/MessageChunks.h"
#include "goby/TooManyHits.h"

using namespace std;

int main (int argc, const char *const argv[]) {
  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  if (argc != 2) {
    cerr << "usage: " << argv[0] << " <basename>" << endl;
    return -1;
  }

  try {
    // this should format the output nicely
    cout.imbue(locale(""));
  } catch (runtime_error& e) {
    // we couldn't set the locale, but that's ok
    cerr << e.what() << endl;
  }
  cout << std::fixed << std::setprecision(6) << boolalpha; 

  const string basename = goby::Alignment::getBasename(argv[1]);
  cout << "Compact Alignment basename = " << basename << endl;

  goby::AlignmentReader alignmentReader = goby::AlignmentReader(basename);
  cout << "Info from header:" << endl;
  cout << "Sorted: " << alignmentReader.isSorted() << endl;
  cout << "Indexed: " << alignmentReader.isIndexed() << endl;

  cout << "Number of target sequences = " << alignmentReader.getNumberOfTargets() << endl;
  
  const vector<unsigned> targetLengths = alignmentReader.getTargetLengths();
  cout << "Number of target length entries = " << targetLengths.size() << endl;
  cout << "smallestSplitQueryIndex = " <<  alignmentReader.getSmallestSplitQueryIndex() << endl;
  cout << "largestSplitQueryIndex = " << alignmentReader.getLargestSplitQueryIndex() << endl;
 
  if (!targetLengths.empty()) {
    cout << "Min target length = " << *min_element(targetLengths.begin(), targetLengths.end()) << endl;
    cout << "Max target length = " << *max_element(targetLengths.begin(), targetLengths.end()) << endl;
    cout << "Mean target length = " << accumulate(targetLengths.begin(), targetLengths.end(), 0) / targetLengths.size() << endl;
  } else {
    cout << "Min target length = 0" << endl;
    cout << "Max target length = 0" << endl;
    cout << "Mean target length = 0" << endl;
  }
  cout << endl;

  cout << "Number of query sequences = " << alignmentReader.getNumberOfQueries() << endl;

  const vector<unsigned> queryLengths = alignmentReader.getQueryLengths();
  cout << "Number of query length entries = " << queryLengths.size() << endl;
  if (!queryLengths.empty()) {
    cout << "Min query length = "  << *min_element(queryLengths.begin(), queryLengths.end()) << endl;
    cout << "Max query length = " << *max_element(queryLengths.begin(), queryLengths.end()) << endl;
    cout << "Mean query length = " << accumulate(queryLengths.begin(), queryLengths.end(), 0) / queryLengths.size() << endl;
  } else {
    cout << "Min query length = 0" << endl;
    cout << "Max query length = 0" << endl;
    cout << "Mean query length = 0" << endl;
  }
  cout << "Constant query lengths = " << alignmentReader.hasConstantQueryLength() << endl;
  cout << "Has query identifiers = " << !alignmentReader.getQueryIdentifiers().empty() << endl;
  cout << "Has target identifiers = " << !alignmentReader.getTargetIdentifiers().empty() << endl;
  cout << endl;

  // Too many hits information
  goby::TooManyHits tmhReader = goby::TooManyHitsReader(basename);
  const vector<unsigned> queryIndicies = tmhReader.getQueryIndicies();
  
  cout << "TMH: aligner threshold = " << tmhReader.getAlignerThreshold() << endl;
  cout << "TMH: number of ambiguous matches = " << queryIndicies.size() << endl;
  cout << "TMH: %ambiguous matches = " << queryIndicies.size() * 100.0f / alignmentReader.getNumberOfQueries() << " %" << endl;

  unsigned maxQueryIndex = 0;
  unsigned maxTargetIndex = 0;
  unsigned numEntries = 0;
  unsigned long numLogicalAlignmentEntries = 0;
  unsigned long total = 0;
  double avgScore = 0;
  unsigned sumNumVariations = 0;

  // from describeAmbigousReads - starts the query index set with the data from the tmh reader
  set<unsigned> alignedQueryIndices(queryIndicies.begin(), queryIndicies.end());

  goby::MessageChunksIterator<goby::AlignmentCollection> alignmentEntriesIterator = alignmentReader.iterator();
  goby::MessageChunksIterator<goby::AlignmentCollection> begin = alignmentEntriesIterator.begin();
  goby::MessageChunksIterator<goby::AlignmentCollection> end = alignmentEntriesIterator.end();
  for (goby::MessageChunksIterator<goby::AlignmentCollection> it = begin; it != end; it++) {
    // cout << "size is " << (*it).alignmententries_size() << endl;
    google::protobuf::RepeatedPtrField<const goby::AlignmentEntry>::const_iterator entryIterator;
    for (entryIterator = (*it).alignmententries().begin(); entryIterator != (*it).alignmententries().end(); entryIterator++) {
      const goby::AlignmentEntry entry = *entryIterator;
      numEntries++;          // Across this file
      numLogicalAlignmentEntries += entry.multiplicity();
      total += entry.query_aligned_length();
      avgScore += entry.score();
      maxQueryIndex = max(maxQueryIndex, entry.query_index());
      maxTargetIndex = max(maxTargetIndex, entry.target_index());
      sumNumVariations += entry.sequence_variations_size();
      alignedQueryIndices.insert(entry.query_index());
    }
  }

  avgScore /= static_cast<double>(numLogicalAlignmentEntries);

  const unsigned numQuerySequences = maxQueryIndex + 1;
  const unsigned numTargetSequences = maxTargetIndex + 1;

  const double avgNumVariationsPerQuery = sumNumVariations / static_cast<double>(numQuerySequences);
        
  cout << "num query indices = " << numQuerySequences << endl;
  cout << "num target indices = " << numTargetSequences << endl;
  cout << "Number of alignment entries = " << numLogicalAlignmentEntries << endl;
  cout << "Number of query indices that matched = " << alignedQueryIndices.size() << endl;
  cout << "Percent matched = " << alignedQueryIndices.size() / static_cast<double>(numQuerySequences) * 100.0f << " %" << endl;

  cout << "Avg query alignment length = " << total / static_cast<double>(numEntries) << endl;
  cout << "Avg score alignment = " << avgScore << endl;
  cout << "Avg number of variations per query sequence = " << avgNumVariationsPerQuery << endl;

  // get the size of the entries file
  const string filename = basename + ".entries";
  struct stat filestatus;
  stat(filename.c_str(), &filestatus);
  const off_t filesize = filestatus.st_size;

  cout << "Average bytes per entry = " << filesize / static_cast<double>(numLogicalAlignmentEntries) << endl;

  // TODO : testing
  /*
  goby::AlignmentWriter alignmentWriter = goby::AlignmentWriter(alignmentReader);
  alignmentWriter.write();

  goby::TooManyHitsWriter tmhWriter = goby::TooManyHitsWriter(tmhReader);
  tmhWriter.write();
  */

  // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
