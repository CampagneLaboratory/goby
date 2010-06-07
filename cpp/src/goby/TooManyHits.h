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

#pragma once

#ifndef GOBY_TOO_MANY_HITS_H
#define GOBY_TOO_MANY_HITS_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Alignments.pb.h"
#include "common.h"

namespace goby {
  // Instantiate classes map<int, int>
  // This does not create an object. It only forces the generation of all
  // of the members of classes vector<int> and vector<char>. It exports
  // them from the DLL and imports them into the .exe file.
  // LIBGOBY_EXPIMP_TEMPLATE template class LIBGOBY_EXPORT std::map<int, int>;

  class LIBGOBY_EXPORT TooManyHits {
    std::string basename;

  protected:
    AlignmentTooManyHits pbTmh;

    // A map from query index to the number of hits at that index.
    std::map<int, int> queryIndex2NumHits;

    // A map from query index to depth/length of match.
    std::map<int, int> queryIndex2Depth;

  public:
    TooManyHits(std::string basename);
    virtual ~TooManyHits(void);

    //TooManyHits(const TooManyHits& from);
    //TooManyHits& operator=(const TooManyHits& from);

    inline const std::string& getBasename() const { return basename; };

    // Number of hits that the aligner considered was too many to report.
    int getAlignerThreshold() const;
    std::vector<int> getQueryIndicies() const;

    int getNumberOfHits(int queryIndex) const;
    int getLengthOfMatch(int queryIndex) const;

    bool isQueryAmbiguous(int queryIndex) const;
    bool isQueryAmbiguous(int queryIndex, int k) const;
    bool isQueryAmbiguous(int queryIndex, int k, int matchLength) const;

    friend std::ostream &operator<<(std::ostream &out, const TooManyHits& tmh);
    friend std::string& operator<<(std::string& out, const TooManyHits& tmh);
    friend TooManyHits& operator<<(TooManyHits& tmh, std::istream& in);
    friend TooManyHits& operator<<(TooManyHits& tmh, std::string& in);
  };

  class LIBGOBY_EXPORT TooManyHitsReader : public TooManyHits {
  public:
    TooManyHitsReader(std::string basename);
    ~TooManyHitsReader(void);
  };

  class LIBGOBY_EXPORT TooManyHitsWriter : public TooManyHits {
  public:
    TooManyHitsWriter(std::string basename);
    ~TooManyHitsWriter(void);

    inline void setAlignerThreshold(int threshold) { pbTmh.set_alignerthreshold(threshold); };

    void write();
  };
}

#endif // GOBY_TOO_MANY_HITS_H
