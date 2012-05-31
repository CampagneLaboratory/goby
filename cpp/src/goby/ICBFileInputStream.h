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

#pragma once

#ifndef ICB_PROTOBUF_GOBY_ICB_FILE_INPUT_STREAM_H__
#define ICB_PROTOBUF_GOBY_ICB_FILE_INPUT_STREAM_H__


#include <string>
#include <iosfwd>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl_lite.h>
#include <google/protobuf/stubs/common.h>

#ifdef _MSC_VER
// Disable Microsoft deprecation warnings for POSIX functions called from this class (open, close)
#pragma warning(push)
#pragma warning(disable:4996)
#endif

using namespace std;

namespace goby {

    // A generic traditional input stream interface.
    //
    // Lots of traditional input streams (e.g. file descriptors, C stdio
    // streams, and C++ iostreams) expose an interface where every read
    // involves copying bytes into a buffer.  If you want to take such an
    // interface and make a ZeroCopyInputStream based on it, simply implement
    // CopyingInputStream and then use CopyingInputStreamAdaptor.
    //
    // CopyingInputStream implementations should avoid buffering if possible.
    // CopyingInputStreamAdaptor does its own buffering and will read data
    // in large blocks.
    class LIBPROTOBUF_EXPORT ICBCopyingInputStream {
        public:
            virtual ~ICBCopyingInputStream();

          // Reads up to "size" bytes into the given buffer.  Returns the number of
          // bytes read.  Read() waits until at least one byte is available, or
          // returns zero if no bytes will ever become available (EOF), or -1 if a
          // permanent read error occurred.
          virtual int Read(void* buffer, int size) = 0;

            // Skips the next "count" bytes of input.  Returns the number of bytes
            // actually skipped.  This will always be exactly equal to "count" unless
            // EOF was reached or a permanent read error occurred.
            //
            // The default implementation just repeatedly calls Read() into a scratch
            // buffer.
            virtual int Skip(int count);

            // Skips the next "count" bytes of input.  Returns the number of bytes
            // actually skipped.  This will always be exactly equal to "count" unless
            // EOF was reached or a permanent read error occurred. Unlike the above
            // Skip() version this takes a std:streamoff which is can handle much bigger
            // files.
            //
            // The default implementation just repeatedly calls Read() into a scratch
            // buffer.
            virtual std::streamoff SkipOffset(std::streamoff count);

    }; // ICBCopyingInputStream

    // ===================================================================

    // A ZeroCopyInputStream which reads from a CopyingInputStream.  This is
    // useful for implementing ZeroCopyInputStreams that read from traditional
    // streams.  Note that this class is not really zero-copy.
    //
    // If you want to read from file descriptors or C++ istreams, this is
    // already implemented for you:  use FileInputStream or IstreamInputStream
    // respectively.
    class LIBPROTOBUF_EXPORT ICBCopyingInputStreamAdaptor : public google::protobuf::io::ZeroCopyInputStream {
        public:
            // Creates a stream that reads from the given CopyingInputStream.
            // If a block_size is given, it specifies the number of bytes that
            // should be read and returned with each call to Next().  Otherwise,
            // a reasonable default is used.  The caller retains ownership of
            // copying_stream unless SetOwnsCopyingStream(true) is called.
            explicit ICBCopyingInputStreamAdaptor(ICBCopyingInputStream* copying_stream,
                                     int block_size = -1);
            ~ICBCopyingInputStreamAdaptor();

            // Call SetOwnsCopyingStream(true) to tell the CopyingInputStreamAdaptor to
            // delete the underlying CopyingInputStream when it is destroyed.
            void SetOwnsCopyingStream(bool value) { owns_copying_stream_ = value; }

            // implements ZeroCopyInputStream ----------------------------------
            bool Next(const void** data, int* size);
            void BackUp(int count);
            bool Skip(int count);
            bool SkipOffset(std::streamoff count);
            google::protobuf::int64 ByteCount() const;

        private:
            // Insures that buffer_ is not NULL.
            void AllocateBufferIfNeeded();
            // Frees the buffer and resets buffer_used_.
            void FreeBuffer();

            // The underlying copying stream.
            ICBCopyingInputStream* copying_stream_;
            bool owns_copying_stream_;

            // True if we have seen a permenant error from the underlying stream.
            bool failed_;

            // The current position of copying_stream_, relative to the point where
            // we started reading.
            google::protobuf::int64 position_;

            // Data is read into this buffer.  It may be NULL if no buffer is currently
            // in use.  Otherwise, it points to an array of size buffer_size_.
            google::protobuf::scoped_array<google::protobuf::uint8> buffer_;
            const int buffer_size_;

            // Number of valid bytes currently in the buffer (i.e. the size last
            // returned by Next()).  0 <= buffer_used_ <= buffer_size_.
            int buffer_used_;

            // Number of bytes in the buffer which were backed up over by a call to
            // BackUp().  These need to be returned again.
            // 0 <= backup_bytes_ <= buffer_used_
            int backup_bytes_;

            GOOGLE_DISALLOW_EVIL_CONSTRUCTORS(ICBCopyingInputStreamAdaptor);
    };

    // ===================================================================

    // A ZeroCopyInputStream which reads from a file descriptor.
    //
    // FileInputStream is preferred over using an ifstream with IstreamInputStream.
    // The latter will introduce an extra layer of buffering, harming performance.
    // Also, it's conceivable that FileInputStream could someday be enhanced
    // to use zero-copy file descriptors on OSs which support them.
    class LIBPROTOBUF_EXPORT ICBFileInputStream : public google::protobuf::io::ZeroCopyInputStream {
        public:
        // Creates a stream that reads from the given Unix file descriptor.
        // If a block_size is given, it specifies the number of bytes that
        // should be read and returned with each call to Next().  Otherwise,
        // a reasonable default is used.
        explicit ICBFileInputStream(int file_descriptor, int block_size = -1);
        ~ICBFileInputStream();

        // Flushes any buffers and closes the underlying file.  Returns false if
        // an error occurs during the process; use GetErrno() to examine the error.
        // Even if an error occurs, the file descriptor is closed when this returns.
        bool Close();

        // By default, the file descriptor is not closed when the stream is
        // destroyed.  Call SetCloseOnDelete(true) to change that.  WARNING:
        // This leaves no way for the caller to detect if close() fails.  If
        // detecting close() errors is important to you, you should arrange
        // to close the descriptor yourself.
        void SetCloseOnDelete(bool value) { copying_input_.SetCloseOnDelete(value); }

        // If an I/O error has occurred on this file descriptor, this is the
        // errno from that error.  Otherwise, this is zero.  Once an error
        // occurs, the stream is broken and all subsequent operations will
        // fail.
        int GetErrno() { return copying_input_.GetErrno(); }

        // implements ZeroCopyInputStream ----------------------------------
        bool Next(const void** data, int* size);
        void BackUp(int count);
        bool Skip(int count);
        bool SkipOffset(std::streamoff count);
        google::protobuf::int64 ByteCount() const;

        private:
            class LIBPROTOBUF_EXPORT ICBCopyingFileInputStream : public ICBCopyingInputStream {
                public:
                    ICBCopyingFileInputStream(int file_descriptor);
                    ~ICBCopyingFileInputStream();

                    bool Close();
                    void SetCloseOnDelete(bool value) { close_on_delete_ = value; }
                    int GetErrno() { return errno_; }

                    // implements CopyingInputStream ---------------------------------
                    int Read(void* buffer, int size);
                    int Skip(int count);
                    std::streamoff SkipOffset(std::streamoff count);

                private:
                    // The file descriptor.
                    const int file_;
                    bool close_on_delete_;
                    bool is_closed_;

                    // The errno of the I/O error, if one has occurred.  Otherwise, zero.
                    int errno_;

                    // Did we try to seek once and fail?  If so, we assume this file descriptor
                    // doesn't support seeking and won't try again.
                    bool previous_seek_failed_;

                    GOOGLE_DISALLOW_EVIL_CONSTRUCTORS(ICBCopyingFileInputStream);
            }; // class ICBCopyingFileInputStream

        ICBCopyingFileInputStream copying_input_;
        ICBCopyingInputStreamAdaptor impl_;

        GOOGLE_DISALLOW_EVIL_CONSTRUCTORS(ICBFileInputStream);
    }; // class ICBFileInputStream

}   // namespace goby



#endif // ICB_PROTOBUF_GOBY_ICB_FILE_INPUT_STREAM_H__
