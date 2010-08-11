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

#ifdef _MSC_VER
#define NOMINMAX // avoid clashing with std::numeric_limits
#endif

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <sys/types.h>
#include <sys/stat.h>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#include "goby/Reads.pb.h"
#include "goby/Reads.h"

using namespace std;

int main(int argc, const char *const argv[]) {
  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  if (argc != 2) {
    cerr << "usage: " << argv[0] << " <basename>" << endl;
    return -1;
  }

  bool verbose = false;
  try {
    // this should format the output nicely
    cout.imbue(locale(""));
  } catch (runtime_error& e) {
    // we couldn't set the locale, but that's ok
    cerr << e.what() << endl;
  }
  cout << std::fixed << std::setprecision(6) << boolalpha; 

  const string basename = goby::Reads::getBasename(argv[1]);
  const string filename = basename + ".compact-reads";
  cout << "Compact reads filename = " << filename << endl;

  // get the size of the file
  struct stat filestatus;
  stat(filename.c_str(), &filestatus);
  const off_t filesize = filestatus.st_size;

  unsigned long long numberOfReads = 0;
  unsigned minReadLength = numeric_limits<unsigned>::max();
  unsigned maxReadLength = numeric_limits<unsigned>::min();
  unsigned long long totalReadLength = 0;
  unsigned long long totalReadLengthPair = 0;

  unsigned long long numberOfIdentifiers = 0;
  unsigned long long numberOfDescriptions = 0;
  unsigned long long numberOfSequences = 0;
  unsigned long long numberOfSequencePairs = 0;
  unsigned long long numberOfQualityScores = 0;
  unsigned long long numberOfQualityScorePairs = 0;

  goby::ReadsReader readsReader = goby::ReadsReader(filename);
  for (goby::ReadEntryIterator it = readsReader.begin(); it != readsReader.end(); it++) {
    goby::ReadEntry entry = *it;
    const unsigned readLength = entry.read_length();

    if (readLength != 0) {
      numberOfReads++;
      totalReadLength += readLength;

      if (entry.has_description()) {
        numberOfDescriptions++;
        if (verbose) {
          cout << "Description found: " << entry.description() << endl;
        }
      }

      if (entry.has_read_identifier()) {
        numberOfIdentifiers++;
        if (verbose) {
          cout << "Identifier found: " << entry.read_identifier() << endl;
        }
      }

      if (entry.has_sequence() && !entry.sequence().empty()) {
        numberOfSequences++;
        if (verbose) {
          cout << "Sequence found: " << entry.sequence() << endl;
        }
      }

      if (entry.has_sequence_pair() && !entry.sequence_pair().empty()) {
        numberOfSequencePairs++;
      }

      if (entry.has_quality_scores() && !entry.quality_scores().empty()) {
        numberOfQualityScores++;
      }

      if (entry.has_quality_scores_pair() && !entry.quality_scores_pair().empty()) {
        numberOfQualityScorePairs++;
      }

      minReadLength = min(minReadLength, readLength);
      maxReadLength = max(maxReadLength, readLength);
    }
  }

  cout << "Average bytes per entry: " << filesize / double(numberOfReads) << endl;
  cout << "Average bytes per base: " << filesize / double(totalReadLength) << endl;
  cout << "Has identifiers = " << (numberOfIdentifiers > 0) << " (" << numberOfIdentifiers << ")" << endl;
  cout << "Has descriptions = " << (numberOfDescriptions > 0) << " (" << numberOfDescriptions << ")" << endl;
  cout << "Has sequences = " << (numberOfSequences > 0) << " (" << numberOfSequences << ")" << endl;
  cout << "Has sequence pairs = " << (numberOfSequencePairs > 0) << " (" << numberOfSequencePairs << ")" << endl;
  cout << "Has quality scores = " << (numberOfQualityScores > 0) << " (" << numberOfQualityScores << ")" << endl;
  cout << "Has quality score pairs = " << (numberOfQualityScorePairs > 0) << " (" << numberOfQualityScorePairs << ")" << endl;
  cout << "Number of entries = " << numberOfReads << endl;
  cout << "Min read length = " << minReadLength << endl;
  cout << "Max read length = " << maxReadLength << endl;
  cout << "Avg read length = " << totalReadLength / (double)numberOfReads << endl;
  cout << "Avg read pair length = " << totalReadLengthPair / (double)numberOfReads << endl;

  // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
