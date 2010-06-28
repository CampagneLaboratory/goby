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
  // Instantiate classes map<unsigned, unsigned>
  // This does not create an object. It only forces the generation of all
  // of the members of classes vector<unsigned> and vector<char>. It exports
  // them from the DLL and imports them into the .exe file.
  // LIBGOBY_EXPIMP_TEMPLATE template class LIBGOBY_EXPORT std::map<unsigned, unsigned>;

  class LIBGOBY_EXPORT TooManyHits {
  protected:
    std::string basename;

    AlignmentTooManyHits pbTmh;

    // A map from query index to the number of hits at that index.
    std::map<unsigned, unsigned> queryIndex2NumHits;

    // A map from query index to depth/length of match.
    std::map<unsigned, unsigned> queryIndex2Depth;

  public:
    TooManyHits(const std::string& basename);
    virtual ~TooManyHits(void);

    inline const std::string& getBasename() const { return basename; };

    // Number of hits that the aligner considered was too many to report.
    unsigned getAlignerThreshold() const;
    std::vector<unsigned> getQueryIndicies() const;

    unsigned getNumberOfHits(unsigned queryIndex) const;
    unsigned getLengthOfMatch(unsigned queryIndex) const;

    bool isQueryAmbiguous(unsigned queryIndex) const;
    bool isQueryAmbiguous(unsigned queryIndex, unsigned k) const;
    bool isQueryAmbiguous(unsigned queryIndex, unsigned k, unsigned matchLength) const;

    friend std::ostream &operator<<(std::ostream &out, const TooManyHits& tmh);
    friend std::string& operator<<(std::string& out, const TooManyHits& tmh);
    friend TooManyHits& operator<<(TooManyHits& tmh, std::istream& in);
    friend TooManyHits& operator<<(TooManyHits& tmh, std::string& in);
  };

  class LIBGOBY_EXPORT TooManyHitsReader : public TooManyHits {
  public:
    TooManyHitsReader(const std::string& basename);
    ~TooManyHitsReader(void);
  };

  class LIBGOBY_EXPORT TooManyHitsWriter : public TooManyHits {
  public:
    TooManyHitsWriter(const std::string& basename);
    TooManyHitsWriter(const std::string& basename, unsigned threshold);
    TooManyHitsWriter(const TooManyHits& tooManyHits);
    ~TooManyHitsWriter(void);

    inline void setAlignerThreshold(unsigned threshold) { pbTmh.set_alignerthreshold(threshold); };

    void append(unsigned queryIndex, unsigned howManyHits, unsigned lengthOfMatch);
    void write();
  };
}

#endif // GOBY_TOO_MANY_HITS_H
