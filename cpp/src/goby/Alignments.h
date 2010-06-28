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

#ifndef GOBY_ALIGNMENTS_H
#define GOBY_ALIGNMENTS_H

#include <string>
#include "common.h"
#include "Alignments.pb.h"
#include "MessageChunks.h"

namespace goby {
  class LIBGOBY_EXPORT Alignment {

  protected:
    std::string basename;
    AlignmentHeader pbHeader;

    // A map of target identifiers (name to index)
    std::map<std::string, unsigned> targetIdentifiers;
    
    // A map of target identifiers (name to index)
    std::map<std::string, unsigned> queryIdentifiers;

  public:
    Alignment(const std::string& basename);
    virtual ~Alignment(void);

    static std::string getBasename(const char* filename);
    static std::string getBasename(const std::string& filename);

    inline const std::string& getBasename() const { return basename; };
    inline unsigned getNumberOfQueries() const { return pbHeader.number_of_queries(); };
    inline unsigned getNumberOfTargets() const { return pbHeader.number_of_targets(); };
    inline unsigned getNumberOfAlignedReads() const { return pbHeader.number_of_aligned_reads(); };
    inline bool hasConstantQueryLength() const { return pbHeader.has_constantquerylength(); };
    inline unsigned getConstantQuerylength() const { return pbHeader.constantquerylength(); };
    inline unsigned getSmallestSplitQueryIndex() const { return pbHeader.smallestsplitqueryindex(); };
    inline unsigned getLargestSplitQueryIndex() const { return  pbHeader.largestsplitqueryindex(); };
    
    std::vector<unsigned> getTargetLengths() const;
    std::vector<unsigned> getQueryLengths() const;   
    
    inline const std::map<std::string, unsigned>& getTargetIdentifiers() const { return targetIdentifiers; };
    inline const std::map<std::string, unsigned>& getQueryIdentifiers() const { return queryIdentifiers; };
  };

  class LIBGOBY_EXPORT AlignmentReader : public Alignment {
    MessageChunksIterator<AlignmentCollection> *messageChunksIterator;
  public:
    AlignmentReader(const std::string& basename);
    ~AlignmentReader(void);
  };

  class LIBGOBY_EXPORT AlignmentWriter : public Alignment {
  public:
    AlignmentWriter(const std::string& basename);
    AlignmentWriter(const Alignment& alignment);
    ~AlignmentWriter(void);

    inline void setNumberOfQueries(unsigned numberOfQueries) { pbHeader.set_number_of_queries(numberOfQueries); };
    inline void setNumberOfTargets(unsigned numberOfTargets) { pbHeader.set_number_of_targets(numberOfTargets); };
    inline void setNumberOfAlignedReads(unsigned numberOfAlignedReads) { pbHeader.set_number_of_aligned_reads(numberOfAlignedReads); };
    inline void setConstantquerylength(unsigned constantQueryLength) { pbHeader.set_constantquerylength(constantQueryLength); };
    inline void setSmallestSplitQueryIndex(unsigned smallestSplitQueryIndex) { pbHeader.set_smallestsplitqueryindex(smallestSplitQueryIndex); };
    inline void setLargestSplitQueryIndex(unsigned largestSplitQueryIndex) { pbHeader.set_largestsplitqueryindex(largestSplitQueryIndex); };

    void setTargetLengths(const std::vector<unsigned>& targetLengths);
    void setTargetLengths(const unsigned* targetLengths);
    // NOTE: Query Lentgh setters are not provided - this information is no longer in the header
    
    // TODO: Target and Query Identifiers
    void write();
  };
}

#endif // GOBY_ALIGNMENTS_H
