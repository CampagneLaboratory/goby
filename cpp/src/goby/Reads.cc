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

#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include "common.h"
#include "MessageChunks.h"
#include "Reads.h"

using namespace std;

namespace goby {
  ReadsIterator::ReadsIterator(const string& filename, const std::streampos position = 0) :
    filename(filename),
    messageChunksIterator(MessageChunksIterator<ReadCollection>(filename, position)),
    readCollection(new ReadCollection),
    currentReadIndex(0) {
  }

  ReadsIterator::ReadsIterator(const ReadsIterator& that) :
    filename(that.filename),
    messageChunksIterator(that.messageChunksIterator),
    readCollection(new ReadCollection),
    currentReadIndex(that.currentReadIndex) {
  }

  ReadsIterator::ReadsIterator(const ReadsIterator& that, const std::streamoff off, const std::ios_base::seekdir dir = std::ios_base::beg) :
    filename(that.filename),
    messageChunksIterator(that.messageChunksIterator, off, dir),
    readCollection(new ReadCollection),
    currentReadIndex(0) {
  }

  ReadsIterator::~ReadsIterator() {
    delete readCollection;
  }

  // Prefix increment operator
  ReadsIterator& ReadsIterator::operator++() {
    ++currentReadIndex;
    // if we're at the end of the current chunk, move on to the next
    if (currentReadIndex >= readCollection->reads_size()) {
      // if there is another chunk, get it otherwise set defaults
      if (messageChunksIterator != messageChunksIterator.end()) {
        messageChunksIterator++;
      } else {
        readCollection->Clear();
      }
      currentReadIndex = 0;
    }
    return *this;
  };

  // Postfix increment operator
  ReadsIterator& ReadsIterator::operator++(int) {
    currentReadIndex++;
    if (currentReadIndex >= readCollection->reads_size()) {
      // if there is another chunk, get it otherwise set defaults
      if (messageChunksIterator != messageChunksIterator.end()) {
        messageChunksIterator++;
      } else {
        readCollection->Clear();
      }
      currentReadIndex = 0;
    }
    return *this;
  };

  bool ReadsIterator::operator==(const ReadsIterator& rhs) const {
    // the filenames must match and the chunk/read indicies must be the same
    return filename == rhs.filename && currentReadIndex == rhs.currentReadIndex && messageChunksIterator == rhs.messageChunksIterator;
  };

  bool ReadsIterator::operator!=(const ReadsIterator& rhs) const {
    // if the filenames or the chunk/read indicies don't match, the reader is different.
    return filename != rhs.filename || currentReadIndex != rhs.currentReadIndex || messageChunksIterator != rhs.messageChunksIterator;
  };

  // return the parsed results for the current chunk
  const ReadEntry& ReadsIterator::operator*() {
    // if we're at the end of the current chunk or at the beginning of a new one
    if (currentReadIndex >= readCollection->reads_size() || currentReadIndex == 0) {
      // if there is another chunk, get it otherwise set defaults
      if (messageChunksIterator != messageChunksIterator.end()) {
        *readCollection = *messageChunksIterator;
      } else {
        readCollection->Clear();
      }
    }
    return readCollection->reads().Get(currentReadIndex);
  };

  ReadEntry* const ReadsIterator::operator->() {
    // TODO
    return NULL;
  };

  ReadsIterator ReadsIterator::begin() const {
    return ReadsIterator(*this);
  };

  ReadsIterator ReadsIterator::end() const {
    return ReadsIterator(*this, static_cast<std::streamoff>(0), std::ios_base::end);
  };


  Reads::Reads(const std::string& basename) : filename(getBasename(basename) + ".compact-reads") {
  }

  Reads::~Reads(void) {
  }

  string Reads::getBasename(const char* filename) {
    return getBasename(string(filename));
  }

  string Reads::getBasename(const string& filename) {
    const string COMPACT_READS_FILE_EXTS[] = {
      ".compact-reads"
    };

    const size_t len = sizeof(COMPACT_READS_FILE_EXTS) / sizeof(COMPACT_READS_FILE_EXTS[0]);
    const size_t dotindex = filename.find_last_of(".");
    if (dotindex != string::npos) {
      const string extension = filename.substr(dotindex);
      for (int i = 0; i < len; i++) {
        if (extension.compare(COMPACT_READS_FILE_EXTS[i]) == 0) {
          return filename.substr(0, dotindex);
        }
      }
    }

    return filename;
  }

  ReadsReader::~ReadsReader(void) {
  }

  ReadsReader::ReadsReader(const string& filename) : Reads(filename) {
  }

  ReadsIterator ReadsReader::iterator() {
    return ReadsIterator(filename);
  }

  ReadsWriter::ReadsWriter(const string& filename, unsigned numberOfEntriesPerChunk) : Reads(getBasename(filename)),
    messageChunksWriter(new MessageChunksWriter<ReadCollection>(filename, numberOfEntriesPerChunk)),
    readCollection(ReadCollection::default_instance()),
    readIndex(0),
    sequence(NULL),
    description(NULL),
    identifier(NULL),
    qualityScores(NULL) {
  }

  ReadsWriter::ReadsWriter(const Reads& reads) : Reads(reads) {
    // TODO: testing only
    this->filename = "foo.compact-reads";
  }

  ReadsWriter::~ReadsWriter(void) {
    readCollection.Clear();
    delete messageChunksWriter;
  }

  void ReadsWriter::appendEntry() {
    // set fields in the new read entry
    goby::ReadEntry *entry = readCollection.add_reads();
    entry->set_read_index(readIndex++);
    if (sequence != NULL) {
      entry->set_sequence(sequence);
      entry->set_read_length(strlen(sequence));
      sequence = NULL;
    } else {
      entry->set_read_length(0);
    }

    if (description != NULL) {
      entry->set_description(description);
      description = NULL;
    }

    if (identifier != NULL) {
      entry->set_read_identifier(identifier);
      identifier = NULL;
    }

    if (qualityScores != NULL) {
      entry->set_quality_scores(qualityScores);
      qualityScores = NULL;
    }

    //cout << entry->DebugString() << endl;

    // and pass it along to the chunk writer
    messageChunksWriter->writeAsNeeded(&readCollection);
  }

  void ReadsWriter::close() {
    messageChunksWriter->flush(&readCollection);
  }
}
