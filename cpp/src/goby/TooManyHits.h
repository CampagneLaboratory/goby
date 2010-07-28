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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>

#include "Alignments.pb.h"
#include "common.h"
#include "hash.h"

namespace goby {
  class LIBGOBY_EXPORT TooManyHits {
  protected:
    std::string basename;

    AlignmentTooManyHits pbTmh;

    // A map from query index to the number of hits at that index.
    LIBGOBY_HASH_MAP<unsigned, unsigned> queryIndex2NumHits;

    // A map from query index to depth/length of match.
    LIBGOBY_HASH_MAP<unsigned, unsigned> queryIndex2Depth;

  public:
    TooManyHits(const std::string& basename);
    virtual ~TooManyHits(void);

    inline const std::string& getBasename() const { return basename; };

    // Number of hits that the aligner considered was too many to report.
    inline unsigned getAlignerThreshold() const { return pbTmh.aligner_threshold(); };
    std::vector<unsigned> getQueryIndicies() const;

    unsigned getNumberOfHits(unsigned queryIndex) const;
    unsigned getLengthOfMatch(unsigned queryIndex) const;

    bool isQueryAmbiguous(unsigned queryIndex) const;
    bool isQueryAmbiguous(unsigned queryIndex, unsigned k) const;
    bool isQueryAmbiguous(unsigned queryIndex, unsigned k, unsigned matchLength) const;
  };

  class LIBGOBY_EXPORT TooManyHitsReader : public TooManyHits {
  public:
    TooManyHitsReader(const std::string& basename);
    ~TooManyHitsReader(void);
  };

  class LIBGOBY_EXPORT TooManyHitsWriter : public TooManyHits {
    bool written;
  public:
    TooManyHitsWriter(const std::string& basename);
    TooManyHitsWriter(const std::string& basename, unsigned threshold);
    ~TooManyHitsWriter(void);

    inline void setAlignerThreshold(unsigned threshold) { pbTmh.set_aligner_threshold(threshold); };

    void append(unsigned queryIndex, unsigned howManyHits, unsigned lengthOfMatch);
    void write();
  };
}

#endif // GOBY_TOO_MANY_HITS_H
