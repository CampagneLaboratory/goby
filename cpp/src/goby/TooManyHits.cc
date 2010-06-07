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
#include <map>
#include <string>
#include <vector>
#include "Alignments.pb.h"
#include "TooManyHits.h"

using namespace std;

namespace goby {
  /*
   * TooManyHits
   */
  TooManyHits::TooManyHits(string basename) {
    // store the basename
    this->basename = basename;
    this->pbTmh = AlignmentTooManyHits::default_instance();
  }

  int TooManyHits::getAlignerThreshold() const {
    return pbTmh.alignerthreshold();
  }
  
  vector<int> TooManyHits::getQueryIndicies() const {
    // copy the keys from the queryIndex2NumHits and return as a new vector 
    vector<int> queryIndicies(queryIndex2NumHits.size());
    int i = 0;
    for (map<int, int>::const_iterator it = queryIndex2NumHits.begin(); it != queryIndex2NumHits.end(); it++) {
      queryIndicies[i++] = it->first;
    }

    return queryIndicies;
  }
  
  int TooManyHits::getNumberOfHits(int queryIndex) const {
    int numberOfHits;
    map<int, int>::const_iterator it = queryIndex2NumHits.find(queryIndex);
    if (it != queryIndex2NumHits.end()) {
      numberOfHits = it->second;
    } else {
      numberOfHits = -1;
    }
    return numberOfHits;
  }

  int TooManyHits::getLengthOfMatch(int queryIndex) const {
    int lengthOfMatch;
    map<int, int>::const_iterator it = queryIndex2Depth.find(queryIndex);
    if (it != queryIndex2NumHits.end()) {
      lengthOfMatch = it->second;
    } else {
      lengthOfMatch = -1;
    }
    return lengthOfMatch;
  }
  
  // Returns true if the query was considered ambiguous by the alignment tool
  bool TooManyHits::isQueryAmbiguous(int queryIndex) const {
    return queryIndex2NumHits.find(queryIndex) != queryIndex2NumHits.end();
  }

  bool TooManyHits::isQueryAmbiguous(int queryIndex, int k) const {
    map<int, int>::const_iterator it = queryIndex2Depth.find(queryIndex);
    if (it == queryIndex2NumHits.end()) {
      return false;
    } if (k >= getAlignerThreshold()) {
      // since k is larger than the aligner threshold, we have to assume the query is
      // ambiguous at k, this is the safe choice.
      return true;
    } else {
      return it->second >= k;
    }
  }

  bool TooManyHits::isQueryAmbiguous(int queryIndex, int k, int matchLength) const {
    if (matchLength < getLengthOfMatch(queryIndex)) {
      return true;
    } else {
      return isQueryAmbiguous(queryIndex, k);
    }
  }

  ostream& operator<<(ostream& out, const TooManyHits& tmh) {
    tmh.pbTmh.SerializeToOstream(&out);
    return out;
  }

  string& operator<<(string& out, const TooManyHits& tmh) {
    tmh.pbTmh.SerializeToString(&out);
    return out;
  }

  TooManyHits& operator<<(TooManyHits& tmh, istream& in) {
    tmh.pbTmh.ParseFromIstream(&in);
    return tmh;
  }

  TooManyHits& operator<<(TooManyHits& tmh, string& in) {
    tmh.pbTmh.ParseFromString(in);
    return tmh;
  }

  TooManyHits::~TooManyHits(void) {
  }

  /*
   * TooManyHitsReader
   */
  TooManyHitsReader::TooManyHitsReader(string basename) : TooManyHits(basename) {
    // open the "tmh" file
    string tmhFilename = basename + ".tmh";
    ifstream tmhStream(tmhFilename.c_str(), ios::in | ios::binary);

    // populate the too many hits object from the file
    if (!pbTmh.ParseFromIstream(&tmhStream)) {
      cerr << "Failed to parse too many hits file" << endl;
    }

    tmhStream.close();

    // populate the query index to number of hits and depth maps
    google::protobuf::RepeatedPtrField<const goby::AmbiguousLocation>::const_iterator hitsIterator;
    for (hitsIterator = pbTmh.hits().begin(); hitsIterator != pbTmh.hits().end(); hitsIterator++) {
      int queryIndex = hitsIterator->query_index();
      queryIndex2NumHits[queryIndex] = hitsIterator->at_least_number_of_hits();
      if (hitsIterator->has_length_of_match()) {
        queryIndex2Depth[queryIndex] = hitsIterator->length_of_match();
      }
    }
  }

  TooManyHitsReader::~TooManyHitsReader(void) {
  }

  /*
   * TooManyHitsWriter
   */
  TooManyHitsWriter::TooManyHitsWriter(string basename) : TooManyHits(basename) {
  }

  TooManyHitsWriter::~TooManyHitsWriter(void) {
  }

}
