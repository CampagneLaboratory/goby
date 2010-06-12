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
#include <vector>

#include "goby/Alignments.h"
#include "goby/TooManyHits.h"

using namespace std;

int main (int argc, const char *const argv[]) {
  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  cout.imbue(std::locale(""));
  cout << std::fixed << std::setprecision(6) << boolalpha; 

  if (argc != 2) {
    cerr << "usage: " << argv[0] << " <basename>" << endl;
    return -1;
  }

  const string basename = goby::Alignment::getBasename(argv[1]);
  cout << "Compact Alignment basename = " << basename << endl;

  goby::Alignment alignmentReader = goby::AlignmentReader(basename);
  cout << "Info from header:" << endl;
  cout << "Number of target sequences = " << alignmentReader.getNumberOfTargets() << endl;
  
  const vector<unsigned> targetLengths = alignmentReader.getTargetLengths();
  cout << "Number of target length entries = " << targetLengths.size() << endl;
  cout << "smallestSplitQueryIndex = " <<  alignmentReader.getSmallestSplitQueryIndex() << endl;
  cout << "largestSplitQueryIndex = " << alignmentReader.getLargestSplitQueryIndex() << endl;
  
  cout << "Min target length = " << *min_element(targetLengths.begin(), targetLengths.end()) << endl;
  cout << "Max target length = " << *max_element(targetLengths.begin(), targetLengths.end()) << endl;
  cout << "Mean target length = TODO" << endl;
  cout << endl;

  cout << "Number of query sequences = " << alignmentReader.getNumberOfQueries() << endl;

  const vector<unsigned> queryLengths = alignmentReader.getQueryLengths();
  cout << "Number of query length entries = " << queryLengths.size() << endl;
  cout << "Min query length = "  << *min_element(queryLengths.begin(), queryLengths.end()) << endl;
  cout << "Max query length = " << *max_element(queryLengths.begin(), queryLengths.end()) << endl;
  cout << "Mean query length = TODO" << endl;
  cout << "Constant query lengths = " << alignmentReader.hasConstantQueryLength() << endl;
  cout << "Has query identifiers = TODO" << endl;
  cout << "Has target identifiers = TODO" << endl;
  cout << endl;

  // Too many hits information
  goby::TooManyHits tmhReader = goby::TooManyHitsReader(basename);
  const vector<unsigned> queryIndicies = tmhReader.getQueryIndicies();
  
  cout << "TMH: aligner threshold = " << tmhReader.getAlignerThreshold() << endl;
  cout << "TMH: number of ambiguous matches = " << queryIndicies.size() << endl;
  cout << "TMH: %ambiguous matches = " << queryIndicies.size() * 100.0f / alignmentReader.getNumberOfQueries() << " %" << endl;

  cout << "num query indices = TODO" << endl;
  cout << "num target indices = TODO" << endl;
  cout << "Number of alignment entries = TODO" << endl;
  cout << "Number of query indices that matched = TODO" << endl;
  cout << "Percent matched = TODO %" << endl;
  cout << "Avg query alignment length = TODO" << endl;
  cout << "Avg score alignment = TODO" << endl;
  cout << "Avg number of variations per query sequence = TODO" << endl;
  cout << "Average bytes per entry = TODO" << endl;

  // TODO : testing
  goby::AlignmentWriter alignmentWriter = goby::AlignmentWriter(alignmentReader);
  alignmentWriter.write();

  goby::TooManyHitsWriter tmhWriter = goby::TooManyHitsWriter(tmhReader);
  tmhWriter.write();

  // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
