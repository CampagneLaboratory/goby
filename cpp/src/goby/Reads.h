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

#ifndef GOBY_READS_H
#define GOBY_READS_H

#include <string>
#include "common.h"
#include "Reads.pb.h"
#include "MessageChunks.h"

namespace goby {
  class LIBGOBY_EXPORT ReadsIterator : public std::iterator<std::input_iterator_tag, ReadEntry> {
    // file name of the goby compact reads file
    std::string filename;

    // iterator over ReadCollections in the compact file
    MessageChunksIterator<ReadCollection> messageChunksIterator;

    // current chunk of read entries
    ReadCollection *readCollection;

    // index of the current read entry in the collection
    int currentReadIndex;

  public:
    ReadsIterator(const std::string& filename, const std::streampos position);
    ReadsIterator(const ReadsIterator& that);
    ReadsIterator(const ReadsIterator& that, std::streamoff off, std::ios_base::seekdir dir);

    virtual ~ReadsIterator();

    // Prefix increment operator
    ReadsIterator& operator++();

    // Postfix increment operator
    ReadsIterator& operator++(int);

    bool operator==(const ReadsIterator& rhs) const;
    bool operator!=(const ReadsIterator& rhs) const;

    const ReadEntry& operator*();
    ReadEntry* const operator->();

    // TODO - remove the operator<< - for testing only
    friend std::ostream &operator<<(std::ostream &out, const ReadsIterator& iter) {
      out << "ostream &operator<< " << iter.readCollection->reads().Get(iter.currentReadIndex).SerializeToOstream(&out);
      return out;
    }

    ReadsIterator begin() const;
    ReadsIterator end() const;
  };

  class LIBGOBY_EXPORT Reads {
  protected:
    std::string filename;

  public:
    Reads(const std::string& basename);
    virtual ~Reads(void);

    static std::string getBasename(const char* filename);
    static std::string getBasename(const std::string& filename);

    inline std::string getBasename() const { return getBasename(filename); };
  };

  class LIBGOBY_EXPORT ReadsReader : public Reads {
  public:
    ReadsReader(const std::string& filename);
    ReadsReader(const ReadsReader& reader);
    ~ReadsReader(void);

    ReadsIterator iterator();
  };

  class LIBGOBY_EXPORT ReadsWriter : public Reads {
    // the underlying message chunk writer
    MessageChunksWriter<ReadCollection> *messageChunksWriter;

    // current chunk of read entries
    ReadCollection readCollection;

    // current read index
    unsigned readIndex;

    char const* sequence;
    char const* description;
    char const* identifier;
    char const* qualityScores;

  public:
    ReadsWriter(const std::string& filename, unsigned numberOfEntriesPerPhunk);
    ReadsWriter(const Reads& reads);
    ~ReadsWriter(void);

    inline void setSequence(char const* sequence) { this->sequence = sequence; };
    inline void setDescription(char const * description) { this->description = description; };
    inline void setIdentifier(char const* identifier) { this->identifier = identifier; };
    inline void setQualityScores(char const* qualityScores) { this->qualityScores = qualityScores; };

    void appendEntry();
    void close();
  };
}

#endif // GOBY_READS_H
