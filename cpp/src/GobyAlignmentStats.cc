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

#include <iostream>
#include <iomanip>
#include <locale>

#include "goby/Alignments.h"
#include "goby/TooManyHits.h"

using namespace std;

int main (int argc, const char *const argv[]) {
  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  cout.imbue(std::locale(""));
  cout << std::fixed << std::setprecision(6); 

  if (argc != 2) {
    cerr << "usage: " << argv[0] << " <basename>" << endl;
    return -1;
  }

  const string basename = goby::Alignments::getBasename(argv[1]);
  cout << "Compact Alignment basename = " << basename << endl;

  goby::TooManyHits tmhReader = goby::TooManyHitsReader(basename);
  const vector<unsigned> queryIndicies = tmhReader.getQueryIndicies();
  
  cout << "TMH: aligner threshold = " << tmhReader.getAlignerThreshold() << endl;
  cout << "TMH: number of ambiguous matches = " << queryIndicies.size() << endl;
  cout << "TMH: %ambiguous matches = " << queryIndicies.size() * 100.0f / 1077455 << " %" << endl;

 // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
