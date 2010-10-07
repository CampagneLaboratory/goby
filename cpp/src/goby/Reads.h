//
// Copyright (C) 2010 Institute for Computational Biomedicine,
//                    Weill Medical College of Cornell University
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#pragma once

#ifndef GOBY_READS_H
#define GOBY_READS_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include "common.h"
#include "Reads.pb.h"
#include "MessageChunks.h"

namespace goby {
  template class LIBGOBY_EXPORT MessageChunksIterator<ReadCollection>;

  class LIBGOBY_EXPORT ReadEntryIterator : public std::iterator<std::input_iterator_tag, ReadEntry> {
    // the file descriptor for the reads file
    int fd;

    // iterator over ReadCollections in the compact file
    MessageChunksIterator<ReadCollection> message_chunks_iterator;

    // the "end" iterator for the collections
    const MessageChunksIterator<ReadCollection> message_chunks_iterator_end;

    // current chunk of read entries
    ReadCollection *read_collection;

    // index of the current read entry in the collection
    int current_read_index;

  public:
    ReadEntryIterator(const int fd, std::streamoff off, std::ios_base::seekdir dir);
    ReadEntryIterator(const ReadEntryIterator& that);

    virtual ~ReadEntryIterator();

    // Prefix increment operator
    ReadEntryIterator& operator++();

    // Postfix increment operator
    ReadEntryIterator& operator++(int);

    bool operator==(const ReadEntryIterator& rhs) const;
    bool operator!=(const ReadEntryIterator& rhs) const;

    const ReadEntry& operator*();
    const ReadEntry* const operator->();
  };

  class LIBGOBY_EXPORT Reads {
  protected:
    // the file name of the reads file
    std::string filename;

  public:
    Reads(const std::string& basename);
    virtual ~Reads(void);

    static std::string getBasename(const char* filename);
    static std::string getBasename(const std::string& filename);

    inline std::string getBasename() const { return getBasename(filename); };
  };

  class LIBGOBY_EXPORT ReadsReader : public Reads {
    // the file descriptor for the reads file
    int fd;

    // the "end" iterator for the read entries
    const ReadEntryIterator *read_entry_iterator_end;

  public:
    ReadsReader(const std::string& filename);
    ReadsReader(const ReadsReader& reader);
    ~ReadsReader(void);

    ReadEntryIterator begin() const;
    ReadEntryIterator* beginPointer() const;
    ReadEntryIterator end() const;
    const ReadEntryIterator* endPointer() const;
  };

  class LIBGOBY_EXPORT ReadsWriter : public Reads {
    // the underlying message chunk writer
    MessageChunksWriter<ReadCollection> *message_chunks_writer;

    // current chunk of read entries
    ReadCollection read_collection;

    // current read index
    unsigned current_read_index;

    char const* sequence;
    char const* description;
    char const* identifier;
    char const* quality_scores;

  public:
    ReadsWriter(const std::string& filename, unsigned number_of_entries_per_chunk = GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK);
    ReadsWriter(const Reads& reads);
    ~ReadsWriter(void);

    inline void setSequence(char const* sequence) { this->sequence = sequence; };
    inline void setDescription(char const * description) { this->description = description; };
    inline void setIdentifier(char const* identifier) { this->identifier = identifier; };
    inline void setQualityScores(char const* quality_scores) { this->quality_scores = quality_scores; };

    void appendEntry();
    void close();
  };
}

#endif // GOBY_READS_H
