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

#include <fstream>
#include <iostream>
#include <string>
#include "Alignments.h"

using namespace std;

namespace goby {
  Alignment::Alignment(string basename) {
    this->basename = basename;
    this->pbHeader = AlignmentHeader::default_instance();
  }

  Alignment::~Alignment(void) {
  }

  string Alignment::getBasename(const char* filename) {
   return getBasename(string(filename));
  }

  string Alignment::getBasename(const string& filename) {
    const string COMPACT_ALIGNMENT_FILE_EXTS[] = {
      ".entries", ".header", ".tmh", ".stats", ".counts"
    };

    const size_t len = sizeof(COMPACT_ALIGNMENT_FILE_EXTS) / sizeof(COMPACT_ALIGNMENT_FILE_EXTS[0]);
    const size_t dotindex = filename.find_last_of(".");
    if (dotindex != string::npos) {
      const string extension = filename.substr(dotindex);
      for (int i = 0; i < len; i++) {
        if (extension.compare(COMPACT_ALIGNMENT_FILE_EXTS[i]) == 0) {
          return filename.substr(0, dotindex);
        }
      }
    }

    return filename;
  }

  AlignmentReader::AlignmentReader(const string& basename) : Alignment(basename) {
    // open the "header" file
    const string headerFilename = basename + ".header";
    ifstream headerStream(headerFilename.c_str(), ios::in | ios::binary);

    // populate the alignment header object from the file
    if (!pbHeader.ParseFromIstream(&headerStream)) {
      cerr << "Failed to parse alignment header file: " << headerFilename << endl;
    }

    headerStream.close();
  }

  AlignmentReader::~AlignmentReader(void) {
  }

}
