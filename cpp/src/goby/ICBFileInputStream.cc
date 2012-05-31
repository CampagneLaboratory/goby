// Protocol Buffers - Google's data interchange format
// Copyright 2008 Google Inc.  All rights reserved.
// http://code.google.com/p/protobuf/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Author: kenton@google.com (Kenton Varda)
// Modifications by: kcd2001@med.cornell.edu (Kevin C. Dorff)
//  Based on original Protocol Buffers design by
//  Sanjay Ghemawat, Jeff Dean, and others.

//
// These classes are based on protobuf's FileInputStream and related classes
// but this implementation adds methods for SkipOffset(std::streamoff)
// to be able to support larger files as the Skip implementation in
// protobuf's versions only support Skip(int).
//

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif
#include <errno.h>
#include <iostream>
#include <algorithm>

#include "ICBFileInputStream.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/stubs/common.h>

using namespace std;

namespace goby {

    namespace {
        // Default block size for Copying{In,Out}putStreamAdaptor.
        static const int kDefaultBlockSize = 8192;
    }  // namespace

    namespace {
        // EINTR sucks.
        int close_no_eintr(int fd) {
          int result;
          do {
            result = close(fd);
          } while (result < 0 && errno == EINTR);
          return result;
        }
    }  // namespace

// ===================================================================

    ICBFileInputStream::ICBFileInputStream(int file_descriptor, int block_size)
            : copying_input_(file_descriptor),
              impl_(&copying_input_, block_size) {
    }

    ICBFileInputStream::~ICBFileInputStream() {}

    bool ICBFileInputStream::Close() {
        return copying_input_.Close();
    }

    bool ICBFileInputStream::Next(const void** data, int* size) {
        return impl_.Next(data, size);
    }

    void ICBFileInputStream::BackUp(int count) {
        impl_.BackUp(count);
    }

    bool ICBFileInputStream::Skip(int count) {
        return impl_.Skip(count);
    }

    bool ICBFileInputStream::SkipOffset(std::streamoff count) {
        return impl_.SkipOffset(count);
    }

    google::protobuf::int64 ICBFileInputStream::ByteCount() const {
        return impl_.ByteCount();
    }

    ICBFileInputStream::ICBCopyingFileInputStream::ICBCopyingFileInputStream(int file_descriptor)
            : file_(file_descriptor),
              close_on_delete_(false),
              is_closed_(false),
              errno_(0),
              previous_seek_failed_(false) {
    }

    ICBFileInputStream::ICBCopyingFileInputStream::~ICBCopyingFileInputStream() {
        if (close_on_delete_) {
            if (!Close()) {
                GOOGLE_LOG(ERROR) << "close() failed: " << strerror(errno_);
            }
        }
    }

    bool ICBFileInputStream::ICBCopyingFileInputStream::Close() {
        GOOGLE_CHECK(!is_closed_);

        is_closed_ = true;
        if (close_no_eintr(file_) != 0) {
            // The docs on close() do not specify whether a file descriptor is still
            // open after close() fails with EIO.  However, the glibc source code
            // seems to indicate that it is not.
            errno_ = errno;
            return false;
        }

        return true;
    }

    int ICBFileInputStream::ICBCopyingFileInputStream::Read(void* buffer, int size) {
        GOOGLE_CHECK(!is_closed_);

        int result;
        do {
            result = read(file_, buffer, size);
        } while (result < 0 && errno == EINTR);

        if (result < 0) {
            // Read error (not EOF).
            errno_ = errno;
        }

        return result;
    }

    int ICBFileInputStream::ICBCopyingFileInputStream::Skip(int count) {
        GOOGLE_CHECK(!is_closed_);

        if (!previous_seek_failed_ &&
            lseek(file_, count, SEEK_CUR) != (off_t)-1) {
            // Seek succeeded.
            return count;
        } else {
            // Failed to seek.

            // Note to self:  Don't seek again.  This file descriptor doesn't
            // support it.
            previous_seek_failed_ = true;

            // Use the default implementation.
            return ICBCopyingInputStream::Skip(count);
        }
    }

    std::streamoff ICBFileInputStream::ICBCopyingFileInputStream::SkipOffset(std::streamoff count) {
        GOOGLE_CHECK(!is_closed_);

        if (!previous_seek_failed_ &&
            lseek(file_, count, SEEK_CUR) != (std::streamoff)-1) {
            // Seek succeeded.
            return count;
        } else {
            // Failed to seek.

            // Note to self:  Don't seek again.  This file descriptor doesn't
            // support it.
            previous_seek_failed_ = true;

            // Use the default implementation.
            return ICBCopyingInputStream::SkipOffset(count);
        }
    }

    // ===================================================================

    ICBCopyingInputStreamAdaptor::ICBCopyingInputStreamAdaptor(ICBCopyingInputStream* copying_stream, int block_size)
            : copying_stream_(copying_stream),
              owns_copying_stream_(false),
              failed_(false),
              position_(0),
              buffer_size_(block_size > 0 ? block_size : kDefaultBlockSize),
              buffer_used_(0),
              backup_bytes_(0) {
    }

    ICBCopyingInputStreamAdaptor::~ICBCopyingInputStreamAdaptor() {
        if (owns_copying_stream_) {
            delete copying_stream_;
        }
    }

    bool ICBCopyingInputStreamAdaptor::Next(const void** data, int* size) {
        if (failed_) {
            // Already failed on a previous read.
            return false;
        }

        AllocateBufferIfNeeded();

        if (backup_bytes_ > 0) {
            // We have data left over from a previous BackUp(), so just return that.
            *data = buffer_.get() + buffer_used_ - backup_bytes_;
            *size = backup_bytes_;
            backup_bytes_ = 0;
            return true;
        }

        // Read new data into the buffer.
        buffer_used_ = copying_stream_->Read(buffer_.get(), buffer_size_);
        if (buffer_used_ <= 0) {
            // EOF or read error.  We don't need the buffer anymore.
            if (buffer_used_ < 0) {
                // Read error (not EOF).
                failed_ = true;
            }
            FreeBuffer();
            return false;
        }
        position_ += buffer_used_;

        *size = buffer_used_;
        *data = buffer_.get();
        return true;
    }

    void ICBCopyingInputStreamAdaptor::BackUp(int count) {
        GOOGLE_CHECK(backup_bytes_ == 0 && buffer_.get() != NULL)
        << " BackUp() can only be called after Next().";
        GOOGLE_CHECK_LE(count, buffer_used_)
        << " Can't back up over more bytes than were returned by the last call"
        " to Next().";
        GOOGLE_CHECK_GE(count, 0)
        << " Parameter to BackUp() can't be negative.";

        backup_bytes_ = count;
    }

    bool ICBCopyingInputStreamAdaptor::Skip(int count) {
        GOOGLE_CHECK_GE(count, 0);

        if (failed_) {
            // Already failed on a previous read.
            return false;
        }

        // First skip any bytes left over from a previous BackUp().
        if (backup_bytes_ >= count) {
            // We have more data left over than we're trying to skip.  Just chop it.
            backup_bytes_ -= count;
            return true;
        }

        count -= backup_bytes_;
        backup_bytes_ = 0;

        int skipped = copying_stream_->Skip(count);
        position_ += skipped;
        return skipped == count;
    }

    bool ICBCopyingInputStreamAdaptor::SkipOffset(std::streamoff count) {
        GOOGLE_CHECK_GE(count, 0);

        if (failed_) {
            // Already failed on a previous read.
            return false;
        }

        // First skip any bytes left over from a previous BackUp().
        if (backup_bytes_ >= count) {
            // We have more data left over than we're trying to skip.  Just chop it.
            backup_bytes_ -= count;
            return true;
        }

        count -= backup_bytes_;
        backup_bytes_ = 0;

        std::streamoff skipped = copying_stream_->SkipOffset(count);
        position_ += skipped;
        return skipped == count;
    }

    google::protobuf::int64 ICBCopyingInputStreamAdaptor::ByteCount() const {
        return position_ - backup_bytes_;
    }

    void ICBCopyingInputStreamAdaptor::AllocateBufferIfNeeded() {
        if (buffer_.get() == NULL) {
            buffer_.reset(new google::protobuf::uint8[buffer_size_]);
        }
    }

    void ICBCopyingInputStreamAdaptor::FreeBuffer() {
        GOOGLE_CHECK_EQ(backup_bytes_, 0);
        buffer_used_ = 0;
        buffer_.reset();
    }

    // ===================================================================

    ICBCopyingInputStream::~ICBCopyingInputStream() {}

    int min(long int x, int y) {
        if (x < y) {
            return x;
        } else {
            return y;
        }
    }

    int ICBCopyingInputStream::Skip(int count) {
      char junk[4096];
      int skipped = 0;
      while (skipped < count) {
        int bytes = Read(junk, min(count - skipped, sizeof(junk)));
        if (bytes <= 0) {
          // EOF or read error.
          return skipped;
        }
        skipped += bytes;
      }
      return skipped;
    }

    std::streamoff ICBCopyingInputStream::SkipOffset(std::streamoff count) {
      char junk[4096];
      std::streamoff skipped = 0;
      while (skipped < count) {
        int bytes = Read(junk, min(count - skipped, sizeof(junk)));
        if (bytes <= 0) {
          // EOF or read error.
          return skipped;
        }
        skipped += bytes;
      }
      return skipped;
    }

} // namespace goby
