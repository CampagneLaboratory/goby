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

#ifndef GOBY_ALIGNMENT_H_
#define GOBY_ALIGNMENT_H_

#include <string>
#include "common.h"
#include "Alignments.pb.h"

namespace goby {
  class LIBGOBY_EXPORT Alignment {

  protected:
    std::string basename;
    AlignmentHeader pbHeader;

  public:
    Alignment(std::string basename);
    virtual ~Alignment(void);

    static std::string getBasename(const char* filename);
    static std::string getBasename(const std::string& filename);

    inline const std::string& getBasename() const { return basename; };
    inline unsigned getNumberOfQueries() const { return pbHeader.number_of_queries(); };
    inline unsigned getNumberOfTargets() const { return pbHeader.number_of_targets(); };
  };

  class LIBGOBY_EXPORT AlignmentReader : public Alignment {
  public:
    AlignmentReader(const std::string& basename);
    ~AlignmentReader(void);
  };

  class LIBGOBY_EXPORT AlignmentWriter : public Alignment {
  public:
    AlignmentWriter(const std::string& basename);
    ~AlignmentWriter(void);

    void write();
  };
}

#endif // _ALIGNMENTS_H_
