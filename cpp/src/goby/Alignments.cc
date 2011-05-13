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

#include <exception>
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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_BOOST_DATE_TIME
#include <boost/date_time/posix_time/posix_time.hpp>
#endif

#ifdef HAVE_BOOST_FILESYSTEM
#include <boost/filesystem/operations.hpp>
#endif

#include "common.h"
#include "hash.h"
#include "Alignments.h"

#ifdef _MSC_VER
// Disable Microsoft deprecation warnings for POSIX functions called from this class (open, close)
#pragma warning(push)
#pragma warning(disable:4996)
#endif

using namespace std;

namespace goby {
  AlignmentEntryIterator::AlignmentEntryIterator(const int fd, streamoff off = 0, ios_base::seekdir dir = ios_base::beg) :
    fd(fd),
    message_chunks_iterator(MessageChunksIterator<AlignmentCollection>(fd, off, 0, dir)),
    message_chunks_iterator_end(MessageChunksIterator<AlignmentCollection>(fd, 0, 0, ios_base::end)),
    alignment_collection(new AlignmentCollection),
    current_alignment_entry_index(0) {
  }

  AlignmentEntryIterator::AlignmentEntryIterator(const AlignmentEntryIterator& that) :
    fd(that.fd),
    message_chunks_iterator(that.message_chunks_iterator),
    message_chunks_iterator_end(that.message_chunks_iterator_end),
    alignment_collection(new AlignmentCollection),
    current_alignment_entry_index(that.current_alignment_entry_index) {
  }

  AlignmentEntryIterator::~AlignmentEntryIterator() {
    delete alignment_collection;
  }

  // Prefix increment operator
  AlignmentEntryIterator& AlignmentEntryIterator::operator++() {
    if (current_alignment_entry_index != -1) {
      ++current_alignment_entry_index;
      // if we're at the end of the current chunk, move on to the next
      if (current_alignment_entry_index >= alignment_collection->alignment_entries_size()) {
        // if there is another chunk, get it otherwise set defaults
        if (message_chunks_iterator != message_chunks_iterator_end) {
          message_chunks_iterator++;
          current_alignment_entry_index = 0;
        } else {
          alignment_collection->Clear();
          current_alignment_entry_index = -1;
        }
      }
    } else {
        cerr << __FILE__ ":" << __LINE__ << " - Attempt to advance past end of fd " << fd << endl;
    }
    return *this;
  };

  // Postfix increment operator
  AlignmentEntryIterator& AlignmentEntryIterator::operator++(int) {
    if (current_alignment_entry_index != -1) {
      current_alignment_entry_index++;
      if (current_alignment_entry_index >= alignment_collection->alignment_entries_size()) {
        // if there is another chunk, get it otherwise set defaults
        if (message_chunks_iterator != message_chunks_iterator_end) {
          message_chunks_iterator++;
          current_alignment_entry_index = 0;
        } else {
          alignment_collection->Clear();
          current_alignment_entry_index = -1;
        }
      }
    } else {
        cerr << __FILE__ ":" << __LINE__ << " - Attempt to advance past end of fd " << fd << endl;
    }
    return *this;
  };

  bool AlignmentEntryIterator::operator==(const AlignmentEntryIterator& rhs) const {
    // the filenames must match and the chunk/read indicies must be the same
    return current_alignment_entry_index == rhs.current_alignment_entry_index && message_chunks_iterator == rhs.message_chunks_iterator;
  };

  bool AlignmentEntryIterator::operator!=(const AlignmentEntryIterator& rhs) const {
    // if the filenames or the chunk/read indicies don't match, the reader is different.
    return current_alignment_entry_index != rhs.current_alignment_entry_index || message_chunks_iterator != rhs.message_chunks_iterator;
  };

  // return the parsed results for the current chunk
  const AlignmentEntry& AlignmentEntryIterator::operator*() {
    // if we're at the end of the current chunk or at the beginning of a new one
    if (current_alignment_entry_index >= alignment_collection->alignment_entries_size() || current_alignment_entry_index == 0) {
      // if there is another chunk, get it otherwise set defaults
      if (message_chunks_iterator != message_chunks_iterator_end) {
        *alignment_collection = *message_chunks_iterator;
      } else {
        alignment_collection->Clear();
      }
    }

    if (alignment_collection->alignment_entries_size() > current_alignment_entry_index) {
      return alignment_collection->alignment_entries().Get(current_alignment_entry_index);
    } else {
      return AlignmentEntry::default_instance();
    }
  };

  const AlignmentEntry* const AlignmentEntryIterator::operator->() {
    // if we're at the end of the current chunk or at the beginning of a new one
    if (current_alignment_entry_index >= alignment_collection->alignment_entries_size() || current_alignment_entry_index == 0) {
      // if there is another chunk, get it otherwise set defaults
      if (message_chunks_iterator != message_chunks_iterator_end) {
        *alignment_collection = *message_chunks_iterator;
      } else {
        alignment_collection->Clear();
      }
    }

    if (alignment_collection->alignment_entries_size() > current_alignment_entry_index) {
      return &alignment_collection->alignment_entries().Get(current_alignment_entry_index);
    } else {
      return &AlignmentEntry::default_instance();
    }
  };

  Alignment::Alignment(const string& basename) : basename(basename),
    header(AlignmentHeader::default_instance()), stats(LIBGOBY_HASH_MAP<string, string>()) {
  }

  Alignment::~Alignment(void) {
  }

  string Alignment::getBasename(const char* filename) {
    return getBasename(string(filename));
  }

  string Alignment::getBasename(const string& filename) {
    const string COMPACT_ALIGNMENT_FILE_EXTS[] = {
      ".entries", ".header", ".tmh", ".stats", ".counts", ".index"
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

  AlignmentReader::AlignmentReader(const string& basename) : Alignment(basename), alignment_entry_iterator_end(NULL) {
    // open the "header" file
    const string header_filename = basename + ".header";
    const int fd = ::open(header_filename.c_str(), O_RDONLY | O_BINARY);

    if (fd > 0) {
      // uncompress file into memory so that it can be parsed
      google::protobuf::io::FileInputStream header_file_stream(fd);
      google::protobuf::io::GzipInputStream gzipHeaderStream(&header_file_stream);

      // populate the alignment header object from the uncompressed data
      if (!header.ParseFromZeroCopyStream(&gzipHeaderStream)) {
        cerr << "Failed to parse alignment header file: " << header_filename << endl;
      }

      // close the streams and files
      header_file_stream.Close();
    } else {
      cerr << "Failed to open alignment header file: " << header_filename << endl;
    }

    // populate the target identifiers
    google::protobuf::RepeatedPtrField<const IdentifierInfo>::const_iterator target_mapping_iterator;
    for (target_mapping_iterator = header.target_name_mapping().mappings().begin(); target_mapping_iterator != header.target_name_mapping().mappings().end(); target_mapping_iterator++) {
      const string target_name = target_mapping_iterator->name();
      const unsigned target_index = target_mapping_iterator->index();
      target_identifiers.insert(pair<string,unsigned>(target_name, target_index));
    }

    // populate the query identifiers
    google::protobuf::RepeatedPtrField<const IdentifierInfo>::const_iterator query_mapping_iterator;
    for (query_mapping_iterator = header.query_name_mapping().mappings().begin(); query_mapping_iterator != header.query_name_mapping().mappings().end(); query_mapping_iterator++) {
      string query_name = query_mapping_iterator->name();
      const unsigned query_index = query_mapping_iterator->index();
      query_identifiers.insert(pair<string,unsigned>(query_name, query_index));
    }

    // populate the query lengths
    if (hasConstantQueryLength()) {
      // The alignment has constant query lengths
      query_lengths.assign(header.number_of_queries(), getConstantQuerylength());
    } else {
      query_lengths.assign(header.query_length().begin(), header.query_length().end());
    }

    // populate the target lengths
    target_lengths.assign(header.target_length().begin(), header.target_length().end());

    // open the "entries" file
    const string entries_filename = basename + ".entries";
    entries_fd = ::open(entries_filename.c_str(), O_RDONLY | O_BINARY);
    if (entries_fd < 0) {
      cerr << "Error opening entries file: " << entries_filename << endl;
    }

    // load the "stats" file
    const string stats_filename = basename + ".stats";
    load_properties(stats_filename, stats);

    alignment_entry_iterator_end = new AlignmentEntryIterator(entries_fd, static_cast<streamoff>(0), ios_base::end);
  }

  AlignmentReader::~AlignmentReader(void) {
    delete alignment_entry_iterator_end;
  }

  AlignmentEntryIterator AlignmentReader::begin() const {
    return AlignmentEntryIterator(entries_fd);
  };

  AlignmentEntryIterator AlignmentReader::end() const {
    return *alignment_entry_iterator_end;
  };

  AlignmentWriter::AlignmentWriter(const string& basename, unsigned number_of_entries_per_chunk) : Alignment(basename),
    entries_chunks_writer(new MessageChunksWriter<AlignmentCollection>(getBasename(basename) + ".entries", number_of_entries_per_chunk)),
    alignment_collection(AlignmentCollection::default_instance()) {
  }

  AlignmentWriter::~AlignmentWriter(void) {
    alignment_collection.Clear();
    delete entries_chunks_writer;

  }
  
  void AlignmentWriter::setTargetLengths(const vector<unsigned>& target_lengths) {
    this->target_lengths = target_lengths;

    // update the header
    header.clear_target_length();
    for (vector<unsigned>::const_iterator it = this->target_lengths.begin(); it < this->target_lengths.end(); it++) {
      header.add_target_length(*it);
    }
  }

  void AlignmentWriter::addTargetLength(const google::protobuf::uint32 targetLength) {
    header.add_target_length(targetLength);
  }

  void AlignmentWriter::setQueryLengthsStoredInEntries(bool value) {
    header.set_query_lengths_stored_in_entries(value);
  }

  void AlignmentWriter::setTargetLengths(const unsigned* target_lengths) {
    int num_elements = sizeof(target_lengths) / sizeof(target_lengths[0]);
    this->target_lengths.resize(num_elements);
    if (num_elements > 0) {
      copy(&target_lengths[0], &target_lengths[num_elements - 1], this->target_lengths.begin());
    }

    // update the header
    header.clear_target_length();
    for (vector<unsigned>::const_iterator it = this->target_lengths.begin(); it < this->target_lengths.end(); it++) {
      header.add_target_length(*it);
    }
  }

  AlignmentEntry* AlignmentWriter::appendEntry() {
    // unlike the reads writer, write any the previous chunks and then return a new entry for the user to populate
    entries_chunks_writer->writeAsNeeded(&alignment_collection);
    return alignment_collection.add_alignment_entries();
  }

  void AlignmentWriter::addTargetIdentifier(const std::string& targetName, const google::protobuf::uint32 targetIndex) {
    vector<google::protobuf::uint32>::iterator found = find(target_name_indexes.begin(), target_name_indexes.end(), targetIndex);
    if (found == target_name_indexes.end()) {
      target_identifiers.insert(pair<string,unsigned>(targetName, targetIndex));
      target_name_indexes.push_back(targetIndex);
      setNumberOfTargets(target_name_indexes.size());
      goby::IdentifierMapping *targetNameMapping = header.mutable_target_name_mapping();
      goby::IdentifierInfo *newMapping = targetNameMapping->add_mappings();
      newMapping->set_name(targetName);
      newMapping->set_index(targetIndex);
    }
  }

  /**
   * If you wish to provide the identifier and have the queryIndex generated
   * automatically. If the identifier has already been registered, you'll get
   * back the same queryIndex as before.
   */
  unsigned AlignmentWriter::addQueryIdentifier(const std::string& queryIdentifier) {
    if (query_identifiers.find(queryIdentifier) == query_identifiers.end()) {
        // New identifier, register it in the local map
        const unsigned newQueryIndex = query_identifiers.size();
        query_identifiers[queryIdentifier] = newQueryIndex;
        // And in the Protobuf map
        goby::IdentifierMapping *queryNameMapping = header.mutable_query_name_mapping();
        goby::IdentifierInfo *newMapping = queryNameMapping->add_mappings();
        newMapping->set_name(queryIdentifier);
        newMapping->set_index(newQueryIndex);
        return newQueryIndex;
    } else {
        // identifier already registered
        return query_identifiers[queryIdentifier];
    }
  }

  /**
   * If you wish to provide the identifier AND the queryIndex. If the 
   * identifier has already been registered, nothing new will happen.
   */
  void AlignmentWriter::addQueryIdentifierWithIndex(const std::string& queryIdentifier, unsigned newQueryIndex) {
    if (query_identifiers.find(queryIdentifier) == query_identifiers.end()) {
        query_identifiers[queryIdentifier] = newQueryIndex;
        // And in the Protobuf map
        goby::IdentifierMapping *queryNameMapping = header.mutable_query_name_mapping();
        goby::IdentifierInfo *newMapping = queryNameMapping->add_mappings();
        newMapping->set_name(queryIdentifier);
        newMapping->set_index(newQueryIndex);
    }
  }

  void AlignmentWriter::close() {
    // Write to the "header" file
    const string headerFilename = basename + ".header";
    int fd = ::open(headerFilename.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, 0644);

    // set up a gzip output stream to compress the header
    google::protobuf::io::FileOutputStream headerFileStream(fd);
    google::protobuf::io::GzipOutputStream gzipHeaderStream(&headerFileStream);

    if (!header.SerializeToZeroCopyStream(&gzipHeaderStream)) {
      cerr << "Failed to write alignment header file: " << headerFilename << endl;
    }

    // close the streams and files
    gzipHeaderStream.Close();
    headerFileStream.Close();    // this call closes the file descriptor as well

    // write the "stats" file
    stats["basename"] = getBasename(basename);
#ifdef HAVE_BOOST_FILESYSTEM
    stats["basename.full"] = boost::filesystem::complete(boost::filesystem::path(basename)).string();
#else
    stats["basename.full"] = basename;
#endif
    stats["min.query.index"] = t_to_string(getSmallestSplitQueryIndex());
    stats["max.query.index"] = t_to_string(getLargestSplitQueryIndex());
    stats["number.of.queries"] = t_to_string(getNumberOfQueries());
    stats["number.alignment.entries"] = t_to_string(getNumberOfAlignedReads());

    const string stats_filename = basename + ".stats";
    ofstream stats_file(stats_filename.c_str(), ios::out | ios::trunc);

    stats_file << "# Statistics for merged alignment." << endl;

#ifdef HAVE_BOOST_DATE_TIME
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    stats_file << "# " << boost::posix_time::to_simple_string(now).c_str() << endl;
#endif

    for (LIBGOBY_HASH_MAP<string, string>::const_iterator it = stats.begin(); it != stats.end(); it++) {
      stats_file << (*it).first << "=" << (*it).second << endl;
    }
    stats_file.close();

    // write any remaining alignment entries
    entries_chunks_writer->close(&alignment_collection);
  }

#ifdef _MSC_VER
#pragma warning(pop)  // Restores the warning state.
#endif
}
