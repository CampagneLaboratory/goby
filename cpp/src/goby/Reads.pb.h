// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: Reads.proto

#ifndef PROTOBUF_Reads_2eproto__INCLUDED
#define PROTOBUF_Reads_2eproto__INCLUDED

#include <string>

#include <google/protobuf/stubs/common.h>

#if GOOGLE_PROTOBUF_VERSION < 2004000
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please update
#error your headers.
#endif
#if 2004001 < GOOGLE_PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/repeated_field.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/generated_message_reflection.h>
// @@protoc_insertion_point(includes)

namespace goby {

// Internal implementation detail -- do not call these.
void  protobuf_AddDesc_Reads_2eproto();
void protobuf_AssignDesc_Reads_2eproto();
void protobuf_ShutdownFile_Reads_2eproto();

class ReadCollection;
class ReadEntry;
class MetaData;

// ===================================================================

class ReadCollection : public ::google::protobuf::Message {
 public:
  ReadCollection();
  virtual ~ReadCollection();
  
  ReadCollection(const ReadCollection& from);
  
  inline ReadCollection& operator=(const ReadCollection& from) {
    CopyFrom(from);
    return *this;
  }
  
  inline const ::google::protobuf::UnknownFieldSet& unknown_fields() const {
    return _unknown_fields_;
  }
  
  inline ::google::protobuf::UnknownFieldSet* mutable_unknown_fields() {
    return &_unknown_fields_;
  }
  
  static const ::google::protobuf::Descriptor* descriptor();
  static const ReadCollection& default_instance();
  
  void Swap(ReadCollection* other);
  
  // implements Message ----------------------------------------------
  
  ReadCollection* New() const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const ReadCollection& from);
  void MergeFrom(const ReadCollection& from);
  void Clear();
  bool IsInitialized() const;
  
  int ByteSize() const;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input);
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const;
  ::google::protobuf::uint8* SerializeWithCachedSizesToArray(::google::protobuf::uint8* output) const;
  int GetCachedSize() const { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const;
  public:
  
  ::google::protobuf::Metadata GetMetadata() const;
  
  // nested types ----------------------------------------------------
  
  // accessors -------------------------------------------------------
  
  // repeated .goby.ReadEntry reads = 1;
  inline int reads_size() const;
  inline void clear_reads();
  static const int kReadsFieldNumber = 1;
  inline const ::goby::ReadEntry& reads(int index) const;
  inline ::goby::ReadEntry* mutable_reads(int index);
  inline ::goby::ReadEntry* add_reads();
  inline const ::google::protobuf::RepeatedPtrField< ::goby::ReadEntry >&
      reads() const;
  inline ::google::protobuf::RepeatedPtrField< ::goby::ReadEntry >*
      mutable_reads();
  
  // @@protoc_insertion_point(class_scope:goby.ReadCollection)
 private:
  
  ::google::protobuf::UnknownFieldSet _unknown_fields_;
  
  ::google::protobuf::RepeatedPtrField< ::goby::ReadEntry > reads_;
  
  mutable int _cached_size_;
  ::google::protobuf::uint32 _has_bits_[(1 + 31) / 32];
  
  friend void  protobuf_AddDesc_Reads_2eproto();
  friend void protobuf_AssignDesc_Reads_2eproto();
  friend void protobuf_ShutdownFile_Reads_2eproto();
  
  void InitAsDefaultInstance();
  static ReadCollection* default_instance_;
};
// -------------------------------------------------------------------

class ReadEntry : public ::google::protobuf::Message {
 public:
  ReadEntry();
  virtual ~ReadEntry();
  
  ReadEntry(const ReadEntry& from);
  
  inline ReadEntry& operator=(const ReadEntry& from) {
    CopyFrom(from);
    return *this;
  }
  
  inline const ::google::protobuf::UnknownFieldSet& unknown_fields() const {
    return _unknown_fields_;
  }
  
  inline ::google::protobuf::UnknownFieldSet* mutable_unknown_fields() {
    return &_unknown_fields_;
  }
  
  static const ::google::protobuf::Descriptor* descriptor();
  static const ReadEntry& default_instance();
  
  void Swap(ReadEntry* other);
  
  // implements Message ----------------------------------------------
  
  ReadEntry* New() const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const ReadEntry& from);
  void MergeFrom(const ReadEntry& from);
  void Clear();
  bool IsInitialized() const;
  
  int ByteSize() const;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input);
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const;
  ::google::protobuf::uint8* SerializeWithCachedSizesToArray(::google::protobuf::uint8* output) const;
  int GetCachedSize() const { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const;
  public:
  
  ::google::protobuf::Metadata GetMetadata() const;
  
  // nested types ----------------------------------------------------
  
  // accessors -------------------------------------------------------
  
  // required uint32 read_index = 1;
  inline bool has_read_index() const;
  inline void clear_read_index();
  static const int kReadIndexFieldNumber = 1;
  inline ::google::protobuf::uint32 read_index() const;
  inline void set_read_index(::google::protobuf::uint32 value);
  
  // optional uint32 barcode_index = 10;
  inline bool has_barcode_index() const;
  inline void clear_barcode_index();
  static const int kBarcodeIndexFieldNumber = 10;
  inline ::google::protobuf::uint32 barcode_index() const;
  inline void set_barcode_index(::google::protobuf::uint32 value);
  
  // optional string read_identifier = 23;
  inline bool has_read_identifier() const;
  inline void clear_read_identifier();
  static const int kReadIdentifierFieldNumber = 23;
  inline const ::std::string& read_identifier() const;
  inline void set_read_identifier(const ::std::string& value);
  inline void set_read_identifier(const char* value);
  inline void set_read_identifier(const char* value, size_t size);
  inline ::std::string* mutable_read_identifier();
  inline ::std::string* release_read_identifier();
  
  // optional string description = 22;
  inline bool has_description() const;
  inline void clear_description();
  static const int kDescriptionFieldNumber = 22;
  inline const ::std::string& description() const;
  inline void set_description(const ::std::string& value);
  inline void set_description(const char* value);
  inline void set_description(const char* value, size_t size);
  inline ::std::string* mutable_description();
  inline ::std::string* release_description();
  
  // required uint32 read_length = 2;
  inline bool has_read_length() const;
  inline void clear_read_length();
  static const int kReadLengthFieldNumber = 2;
  inline ::google::protobuf::uint32 read_length() const;
  inline void set_read_length(::google::protobuf::uint32 value);
  
  // optional bytes sequence = 3;
  inline bool has_sequence() const;
  inline void clear_sequence();
  static const int kSequenceFieldNumber = 3;
  inline const ::std::string& sequence() const;
  inline void set_sequence(const ::std::string& value);
  inline void set_sequence(const char* value);
  inline void set_sequence(const void* value, size_t size);
  inline ::std::string* mutable_sequence();
  inline ::std::string* release_sequence();
  
  // optional bytes sequence_pair = 5;
  inline bool has_sequence_pair() const;
  inline void clear_sequence_pair();
  static const int kSequencePairFieldNumber = 5;
  inline const ::std::string& sequence_pair() const;
  inline void set_sequence_pair(const ::std::string& value);
  inline void set_sequence_pair(const char* value);
  inline void set_sequence_pair(const void* value, size_t size);
  inline ::std::string* mutable_sequence_pair();
  inline ::std::string* release_sequence_pair();
  
  // optional uint32 read_length_pair = 6;
  inline bool has_read_length_pair() const;
  inline void clear_read_length_pair();
  static const int kReadLengthPairFieldNumber = 6;
  inline ::google::protobuf::uint32 read_length_pair() const;
  inline void set_read_length_pair(::google::protobuf::uint32 value);
  
  // optional bytes quality_scores = 4;
  inline bool has_quality_scores() const;
  inline void clear_quality_scores();
  static const int kQualityScoresFieldNumber = 4;
  inline const ::std::string& quality_scores() const;
  inline void set_quality_scores(const ::std::string& value);
  inline void set_quality_scores(const char* value);
  inline void set_quality_scores(const void* value, size_t size);
  inline ::std::string* mutable_quality_scores();
  inline ::std::string* release_quality_scores();
  
  // optional bytes quality_scores_pair = 7;
  inline bool has_quality_scores_pair() const;
  inline void clear_quality_scores_pair();
  static const int kQualityScoresPairFieldNumber = 7;
  inline const ::std::string& quality_scores_pair() const;
  inline void set_quality_scores_pair(const ::std::string& value);
  inline void set_quality_scores_pair(const char* value);
  inline void set_quality_scores_pair(const void* value, size_t size);
  inline ::std::string* mutable_quality_scores_pair();
  inline ::std::string* release_quality_scores_pair();
  
  // optional bytes compressed_data = 8;
  inline bool has_compressed_data() const;
  inline void clear_compressed_data();
  static const int kCompressedDataFieldNumber = 8;
  inline const ::std::string& compressed_data() const;
  inline void set_compressed_data(const ::std::string& value);
  inline void set_compressed_data(const char* value);
  inline void set_compressed_data(const void* value, size_t size);
  inline ::std::string* mutable_compressed_data();
  inline ::std::string* release_compressed_data();
  
  // repeated .goby.MetaData meta_data = 25;
  inline int meta_data_size() const;
  inline void clear_meta_data();
  static const int kMetaDataFieldNumber = 25;
  inline const ::goby::MetaData& meta_data(int index) const;
  inline ::goby::MetaData* mutable_meta_data(int index);
  inline ::goby::MetaData* add_meta_data();
  inline const ::google::protobuf::RepeatedPtrField< ::goby::MetaData >&
      meta_data() const;
  inline ::google::protobuf::RepeatedPtrField< ::goby::MetaData >*
      mutable_meta_data();
  
  // @@protoc_insertion_point(class_scope:goby.ReadEntry)
 private:
  inline void set_has_read_index();
  inline void clear_has_read_index();
  inline void set_has_barcode_index();
  inline void clear_has_barcode_index();
  inline void set_has_read_identifier();
  inline void clear_has_read_identifier();
  inline void set_has_description();
  inline void clear_has_description();
  inline void set_has_read_length();
  inline void clear_has_read_length();
  inline void set_has_sequence();
  inline void clear_has_sequence();
  inline void set_has_sequence_pair();
  inline void clear_has_sequence_pair();
  inline void set_has_read_length_pair();
  inline void clear_has_read_length_pair();
  inline void set_has_quality_scores();
  inline void clear_has_quality_scores();
  inline void set_has_quality_scores_pair();
  inline void clear_has_quality_scores_pair();
  inline void set_has_compressed_data();
  inline void clear_has_compressed_data();
  
  ::google::protobuf::UnknownFieldSet _unknown_fields_;
  
  ::google::protobuf::uint32 read_index_;
  ::google::protobuf::uint32 barcode_index_;
  ::std::string* read_identifier_;
  ::std::string* description_;
  ::std::string* sequence_;
  ::google::protobuf::uint32 read_length_;
  ::google::protobuf::uint32 read_length_pair_;
  ::std::string* sequence_pair_;
  ::std::string* quality_scores_;
  ::std::string* quality_scores_pair_;
  ::std::string* compressed_data_;
  ::google::protobuf::RepeatedPtrField< ::goby::MetaData > meta_data_;
  
  mutable int _cached_size_;
  ::google::protobuf::uint32 _has_bits_[(12 + 31) / 32];
  
  friend void  protobuf_AddDesc_Reads_2eproto();
  friend void protobuf_AssignDesc_Reads_2eproto();
  friend void protobuf_ShutdownFile_Reads_2eproto();
  
  void InitAsDefaultInstance();
  static ReadEntry* default_instance_;
};
// -------------------------------------------------------------------

class MetaData : public ::google::protobuf::Message {
 public:
  MetaData();
  virtual ~MetaData();
  
  MetaData(const MetaData& from);
  
  inline MetaData& operator=(const MetaData& from) {
    CopyFrom(from);
    return *this;
  }
  
  inline const ::google::protobuf::UnknownFieldSet& unknown_fields() const {
    return _unknown_fields_;
  }
  
  inline ::google::protobuf::UnknownFieldSet* mutable_unknown_fields() {
    return &_unknown_fields_;
  }
  
  static const ::google::protobuf::Descriptor* descriptor();
  static const MetaData& default_instance();
  
  void Swap(MetaData* other);
  
  // implements Message ----------------------------------------------
  
  MetaData* New() const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const MetaData& from);
  void MergeFrom(const MetaData& from);
  void Clear();
  bool IsInitialized() const;
  
  int ByteSize() const;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input);
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const;
  ::google::protobuf::uint8* SerializeWithCachedSizesToArray(::google::protobuf::uint8* output) const;
  int GetCachedSize() const { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const;
  public:
  
  ::google::protobuf::Metadata GetMetadata() const;
  
  // nested types ----------------------------------------------------
  
  // accessors -------------------------------------------------------
  
  // required string key = 1;
  inline bool has_key() const;
  inline void clear_key();
  static const int kKeyFieldNumber = 1;
  inline const ::std::string& key() const;
  inline void set_key(const ::std::string& value);
  inline void set_key(const char* value);
  inline void set_key(const char* value, size_t size);
  inline ::std::string* mutable_key();
  inline ::std::string* release_key();
  
  // required string value = 2;
  inline bool has_value() const;
  inline void clear_value();
  static const int kValueFieldNumber = 2;
  inline const ::std::string& value() const;
  inline void set_value(const ::std::string& value);
  inline void set_value(const char* value);
  inline void set_value(const char* value, size_t size);
  inline ::std::string* mutable_value();
  inline ::std::string* release_value();
  
  // @@protoc_insertion_point(class_scope:goby.MetaData)
 private:
  inline void set_has_key();
  inline void clear_has_key();
  inline void set_has_value();
  inline void clear_has_value();
  
  ::google::protobuf::UnknownFieldSet _unknown_fields_;
  
  ::std::string* key_;
  ::std::string* value_;
  
  mutable int _cached_size_;
  ::google::protobuf::uint32 _has_bits_[(2 + 31) / 32];
  
  friend void  protobuf_AddDesc_Reads_2eproto();
  friend void protobuf_AssignDesc_Reads_2eproto();
  friend void protobuf_ShutdownFile_Reads_2eproto();
  
  void InitAsDefaultInstance();
  static MetaData* default_instance_;
};
// ===================================================================


// ===================================================================

// ReadCollection

// repeated .goby.ReadEntry reads = 1;
inline int ReadCollection::reads_size() const {
  return reads_.size();
}
inline void ReadCollection::clear_reads() {
  reads_.Clear();
}
inline const ::goby::ReadEntry& ReadCollection::reads(int index) const {
  return reads_.Get(index);
}
inline ::goby::ReadEntry* ReadCollection::mutable_reads(int index) {
  return reads_.Mutable(index);
}
inline ::goby::ReadEntry* ReadCollection::add_reads() {
  return reads_.Add();
}
inline const ::google::protobuf::RepeatedPtrField< ::goby::ReadEntry >&
ReadCollection::reads() const {
  return reads_;
}
inline ::google::protobuf::RepeatedPtrField< ::goby::ReadEntry >*
ReadCollection::mutable_reads() {
  return &reads_;
}

// -------------------------------------------------------------------

// ReadEntry

// required uint32 read_index = 1;
inline bool ReadEntry::has_read_index() const {
  return (_has_bits_[0] & 0x00000001u) != 0;
}
inline void ReadEntry::set_has_read_index() {
  _has_bits_[0] |= 0x00000001u;
}
inline void ReadEntry::clear_has_read_index() {
  _has_bits_[0] &= ~0x00000001u;
}
inline void ReadEntry::clear_read_index() {
  read_index_ = 0u;
  clear_has_read_index();
}
inline ::google::protobuf::uint32 ReadEntry::read_index() const {
  return read_index_;
}
inline void ReadEntry::set_read_index(::google::protobuf::uint32 value) {
  set_has_read_index();
  read_index_ = value;
}

// optional uint32 barcode_index = 10;
inline bool ReadEntry::has_barcode_index() const {
  return (_has_bits_[0] & 0x00000002u) != 0;
}
inline void ReadEntry::set_has_barcode_index() {
  _has_bits_[0] |= 0x00000002u;
}
inline void ReadEntry::clear_has_barcode_index() {
  _has_bits_[0] &= ~0x00000002u;
}
inline void ReadEntry::clear_barcode_index() {
  barcode_index_ = 0u;
  clear_has_barcode_index();
}
inline ::google::protobuf::uint32 ReadEntry::barcode_index() const {
  return barcode_index_;
}
inline void ReadEntry::set_barcode_index(::google::protobuf::uint32 value) {
  set_has_barcode_index();
  barcode_index_ = value;
}

// optional string read_identifier = 23;
inline bool ReadEntry::has_read_identifier() const {
  return (_has_bits_[0] & 0x00000004u) != 0;
}
inline void ReadEntry::set_has_read_identifier() {
  _has_bits_[0] |= 0x00000004u;
}
inline void ReadEntry::clear_has_read_identifier() {
  _has_bits_[0] &= ~0x00000004u;
}
inline void ReadEntry::clear_read_identifier() {
  if (read_identifier_ != &::google::protobuf::internal::kEmptyString) {
    read_identifier_->clear();
  }
  clear_has_read_identifier();
}
inline const ::std::string& ReadEntry::read_identifier() const {
  return *read_identifier_;
}
inline void ReadEntry::set_read_identifier(const ::std::string& value) {
  set_has_read_identifier();
  if (read_identifier_ == &::google::protobuf::internal::kEmptyString) {
    read_identifier_ = new ::std::string;
  }
  read_identifier_->assign(value);
}
inline void ReadEntry::set_read_identifier(const char* value) {
  set_has_read_identifier();
  if (read_identifier_ == &::google::protobuf::internal::kEmptyString) {
    read_identifier_ = new ::std::string;
  }
  read_identifier_->assign(value);
}
inline void ReadEntry::set_read_identifier(const char* value, size_t size) {
  set_has_read_identifier();
  if (read_identifier_ == &::google::protobuf::internal::kEmptyString) {
    read_identifier_ = new ::std::string;
  }
  read_identifier_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* ReadEntry::mutable_read_identifier() {
  set_has_read_identifier();
  if (read_identifier_ == &::google::protobuf::internal::kEmptyString) {
    read_identifier_ = new ::std::string;
  }
  return read_identifier_;
}
inline ::std::string* ReadEntry::release_read_identifier() {
  clear_has_read_identifier();
  if (read_identifier_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = read_identifier_;
    read_identifier_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}

// optional string description = 22;
inline bool ReadEntry::has_description() const {
  return (_has_bits_[0] & 0x00000008u) != 0;
}
inline void ReadEntry::set_has_description() {
  _has_bits_[0] |= 0x00000008u;
}
inline void ReadEntry::clear_has_description() {
  _has_bits_[0] &= ~0x00000008u;
}
inline void ReadEntry::clear_description() {
  if (description_ != &::google::protobuf::internal::kEmptyString) {
    description_->clear();
  }
  clear_has_description();
}
inline const ::std::string& ReadEntry::description() const {
  return *description_;
}
inline void ReadEntry::set_description(const ::std::string& value) {
  set_has_description();
  if (description_ == &::google::protobuf::internal::kEmptyString) {
    description_ = new ::std::string;
  }
  description_->assign(value);
}
inline void ReadEntry::set_description(const char* value) {
  set_has_description();
  if (description_ == &::google::protobuf::internal::kEmptyString) {
    description_ = new ::std::string;
  }
  description_->assign(value);
}
inline void ReadEntry::set_description(const char* value, size_t size) {
  set_has_description();
  if (description_ == &::google::protobuf::internal::kEmptyString) {
    description_ = new ::std::string;
  }
  description_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* ReadEntry::mutable_description() {
  set_has_description();
  if (description_ == &::google::protobuf::internal::kEmptyString) {
    description_ = new ::std::string;
  }
  return description_;
}
inline ::std::string* ReadEntry::release_description() {
  clear_has_description();
  if (description_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = description_;
    description_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}

// required uint32 read_length = 2;
inline bool ReadEntry::has_read_length() const {
  return (_has_bits_[0] & 0x00000010u) != 0;
}
inline void ReadEntry::set_has_read_length() {
  _has_bits_[0] |= 0x00000010u;
}
inline void ReadEntry::clear_has_read_length() {
  _has_bits_[0] &= ~0x00000010u;
}
inline void ReadEntry::clear_read_length() {
  read_length_ = 0u;
  clear_has_read_length();
}
inline ::google::protobuf::uint32 ReadEntry::read_length() const {
  return read_length_;
}
inline void ReadEntry::set_read_length(::google::protobuf::uint32 value) {
  set_has_read_length();
  read_length_ = value;
}

// optional bytes sequence = 3;
inline bool ReadEntry::has_sequence() const {
  return (_has_bits_[0] & 0x00000020u) != 0;
}
inline void ReadEntry::set_has_sequence() {
  _has_bits_[0] |= 0x00000020u;
}
inline void ReadEntry::clear_has_sequence() {
  _has_bits_[0] &= ~0x00000020u;
}
inline void ReadEntry::clear_sequence() {
  if (sequence_ != &::google::protobuf::internal::kEmptyString) {
    sequence_->clear();
  }
  clear_has_sequence();
}
inline const ::std::string& ReadEntry::sequence() const {
  return *sequence_;
}
inline void ReadEntry::set_sequence(const ::std::string& value) {
  set_has_sequence();
  if (sequence_ == &::google::protobuf::internal::kEmptyString) {
    sequence_ = new ::std::string;
  }
  sequence_->assign(value);
}
inline void ReadEntry::set_sequence(const char* value) {
  set_has_sequence();
  if (sequence_ == &::google::protobuf::internal::kEmptyString) {
    sequence_ = new ::std::string;
  }
  sequence_->assign(value);
}
inline void ReadEntry::set_sequence(const void* value, size_t size) {
  set_has_sequence();
  if (sequence_ == &::google::protobuf::internal::kEmptyString) {
    sequence_ = new ::std::string;
  }
  sequence_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* ReadEntry::mutable_sequence() {
  set_has_sequence();
  if (sequence_ == &::google::protobuf::internal::kEmptyString) {
    sequence_ = new ::std::string;
  }
  return sequence_;
}
inline ::std::string* ReadEntry::release_sequence() {
  clear_has_sequence();
  if (sequence_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = sequence_;
    sequence_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}

// optional bytes sequence_pair = 5;
inline bool ReadEntry::has_sequence_pair() const {
  return (_has_bits_[0] & 0x00000040u) != 0;
}
inline void ReadEntry::set_has_sequence_pair() {
  _has_bits_[0] |= 0x00000040u;
}
inline void ReadEntry::clear_has_sequence_pair() {
  _has_bits_[0] &= ~0x00000040u;
}
inline void ReadEntry::clear_sequence_pair() {
  if (sequence_pair_ != &::google::protobuf::internal::kEmptyString) {
    sequence_pair_->clear();
  }
  clear_has_sequence_pair();
}
inline const ::std::string& ReadEntry::sequence_pair() const {
  return *sequence_pair_;
}
inline void ReadEntry::set_sequence_pair(const ::std::string& value) {
  set_has_sequence_pair();
  if (sequence_pair_ == &::google::protobuf::internal::kEmptyString) {
    sequence_pair_ = new ::std::string;
  }
  sequence_pair_->assign(value);
}
inline void ReadEntry::set_sequence_pair(const char* value) {
  set_has_sequence_pair();
  if (sequence_pair_ == &::google::protobuf::internal::kEmptyString) {
    sequence_pair_ = new ::std::string;
  }
  sequence_pair_->assign(value);
}
inline void ReadEntry::set_sequence_pair(const void* value, size_t size) {
  set_has_sequence_pair();
  if (sequence_pair_ == &::google::protobuf::internal::kEmptyString) {
    sequence_pair_ = new ::std::string;
  }
  sequence_pair_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* ReadEntry::mutable_sequence_pair() {
  set_has_sequence_pair();
  if (sequence_pair_ == &::google::protobuf::internal::kEmptyString) {
    sequence_pair_ = new ::std::string;
  }
  return sequence_pair_;
}
inline ::std::string* ReadEntry::release_sequence_pair() {
  clear_has_sequence_pair();
  if (sequence_pair_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = sequence_pair_;
    sequence_pair_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}

// optional uint32 read_length_pair = 6;
inline bool ReadEntry::has_read_length_pair() const {
  return (_has_bits_[0] & 0x00000080u) != 0;
}
inline void ReadEntry::set_has_read_length_pair() {
  _has_bits_[0] |= 0x00000080u;
}
inline void ReadEntry::clear_has_read_length_pair() {
  _has_bits_[0] &= ~0x00000080u;
}
inline void ReadEntry::clear_read_length_pair() {
  read_length_pair_ = 0u;
  clear_has_read_length_pair();
}
inline ::google::protobuf::uint32 ReadEntry::read_length_pair() const {
  return read_length_pair_;
}
inline void ReadEntry::set_read_length_pair(::google::protobuf::uint32 value) {
  set_has_read_length_pair();
  read_length_pair_ = value;
}

// optional bytes quality_scores = 4;
inline bool ReadEntry::has_quality_scores() const {
  return (_has_bits_[0] & 0x00000100u) != 0;
}
inline void ReadEntry::set_has_quality_scores() {
  _has_bits_[0] |= 0x00000100u;
}
inline void ReadEntry::clear_has_quality_scores() {
  _has_bits_[0] &= ~0x00000100u;
}
inline void ReadEntry::clear_quality_scores() {
  if (quality_scores_ != &::google::protobuf::internal::kEmptyString) {
    quality_scores_->clear();
  }
  clear_has_quality_scores();
}
inline const ::std::string& ReadEntry::quality_scores() const {
  return *quality_scores_;
}
inline void ReadEntry::set_quality_scores(const ::std::string& value) {
  set_has_quality_scores();
  if (quality_scores_ == &::google::protobuf::internal::kEmptyString) {
    quality_scores_ = new ::std::string;
  }
  quality_scores_->assign(value);
}
inline void ReadEntry::set_quality_scores(const char* value) {
  set_has_quality_scores();
  if (quality_scores_ == &::google::protobuf::internal::kEmptyString) {
    quality_scores_ = new ::std::string;
  }
  quality_scores_->assign(value);
}
inline void ReadEntry::set_quality_scores(const void* value, size_t size) {
  set_has_quality_scores();
  if (quality_scores_ == &::google::protobuf::internal::kEmptyString) {
    quality_scores_ = new ::std::string;
  }
  quality_scores_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* ReadEntry::mutable_quality_scores() {
  set_has_quality_scores();
  if (quality_scores_ == &::google::protobuf::internal::kEmptyString) {
    quality_scores_ = new ::std::string;
  }
  return quality_scores_;
}
inline ::std::string* ReadEntry::release_quality_scores() {
  clear_has_quality_scores();
  if (quality_scores_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = quality_scores_;
    quality_scores_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}

// optional bytes quality_scores_pair = 7;
inline bool ReadEntry::has_quality_scores_pair() const {
  return (_has_bits_[0] & 0x00000200u) != 0;
}
inline void ReadEntry::set_has_quality_scores_pair() {
  _has_bits_[0] |= 0x00000200u;
}
inline void ReadEntry::clear_has_quality_scores_pair() {
  _has_bits_[0] &= ~0x00000200u;
}
inline void ReadEntry::clear_quality_scores_pair() {
  if (quality_scores_pair_ != &::google::protobuf::internal::kEmptyString) {
    quality_scores_pair_->clear();
  }
  clear_has_quality_scores_pair();
}
inline const ::std::string& ReadEntry::quality_scores_pair() const {
  return *quality_scores_pair_;
}
inline void ReadEntry::set_quality_scores_pair(const ::std::string& value) {
  set_has_quality_scores_pair();
  if (quality_scores_pair_ == &::google::protobuf::internal::kEmptyString) {
    quality_scores_pair_ = new ::std::string;
  }
  quality_scores_pair_->assign(value);
}
inline void ReadEntry::set_quality_scores_pair(const char* value) {
  set_has_quality_scores_pair();
  if (quality_scores_pair_ == &::google::protobuf::internal::kEmptyString) {
    quality_scores_pair_ = new ::std::string;
  }
  quality_scores_pair_->assign(value);
}
inline void ReadEntry::set_quality_scores_pair(const void* value, size_t size) {
  set_has_quality_scores_pair();
  if (quality_scores_pair_ == &::google::protobuf::internal::kEmptyString) {
    quality_scores_pair_ = new ::std::string;
  }
  quality_scores_pair_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* ReadEntry::mutable_quality_scores_pair() {
  set_has_quality_scores_pair();
  if (quality_scores_pair_ == &::google::protobuf::internal::kEmptyString) {
    quality_scores_pair_ = new ::std::string;
  }
  return quality_scores_pair_;
}
inline ::std::string* ReadEntry::release_quality_scores_pair() {
  clear_has_quality_scores_pair();
  if (quality_scores_pair_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = quality_scores_pair_;
    quality_scores_pair_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}

// optional bytes compressed_data = 8;
inline bool ReadEntry::has_compressed_data() const {
  return (_has_bits_[0] & 0x00000400u) != 0;
}
inline void ReadEntry::set_has_compressed_data() {
  _has_bits_[0] |= 0x00000400u;
}
inline void ReadEntry::clear_has_compressed_data() {
  _has_bits_[0] &= ~0x00000400u;
}
inline void ReadEntry::clear_compressed_data() {
  if (compressed_data_ != &::google::protobuf::internal::kEmptyString) {
    compressed_data_->clear();
  }
  clear_has_compressed_data();
}
inline const ::std::string& ReadEntry::compressed_data() const {
  return *compressed_data_;
}
inline void ReadEntry::set_compressed_data(const ::std::string& value) {
  set_has_compressed_data();
  if (compressed_data_ == &::google::protobuf::internal::kEmptyString) {
    compressed_data_ = new ::std::string;
  }
  compressed_data_->assign(value);
}
inline void ReadEntry::set_compressed_data(const char* value) {
  set_has_compressed_data();
  if (compressed_data_ == &::google::protobuf::internal::kEmptyString) {
    compressed_data_ = new ::std::string;
  }
  compressed_data_->assign(value);
}
inline void ReadEntry::set_compressed_data(const void* value, size_t size) {
  set_has_compressed_data();
  if (compressed_data_ == &::google::protobuf::internal::kEmptyString) {
    compressed_data_ = new ::std::string;
  }
  compressed_data_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* ReadEntry::mutable_compressed_data() {
  set_has_compressed_data();
  if (compressed_data_ == &::google::protobuf::internal::kEmptyString) {
    compressed_data_ = new ::std::string;
  }
  return compressed_data_;
}
inline ::std::string* ReadEntry::release_compressed_data() {
  clear_has_compressed_data();
  if (compressed_data_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = compressed_data_;
    compressed_data_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}

// repeated .goby.MetaData meta_data = 25;
inline int ReadEntry::meta_data_size() const {
  return meta_data_.size();
}
inline void ReadEntry::clear_meta_data() {
  meta_data_.Clear();
}
inline const ::goby::MetaData& ReadEntry::meta_data(int index) const {
  return meta_data_.Get(index);
}
inline ::goby::MetaData* ReadEntry::mutable_meta_data(int index) {
  return meta_data_.Mutable(index);
}
inline ::goby::MetaData* ReadEntry::add_meta_data() {
  return meta_data_.Add();
}
inline const ::google::protobuf::RepeatedPtrField< ::goby::MetaData >&
ReadEntry::meta_data() const {
  return meta_data_;
}
inline ::google::protobuf::RepeatedPtrField< ::goby::MetaData >*
ReadEntry::mutable_meta_data() {
  return &meta_data_;
}

// -------------------------------------------------------------------

// MetaData

// required string key = 1;
inline bool MetaData::has_key() const {
  return (_has_bits_[0] & 0x00000001u) != 0;
}
inline void MetaData::set_has_key() {
  _has_bits_[0] |= 0x00000001u;
}
inline void MetaData::clear_has_key() {
  _has_bits_[0] &= ~0x00000001u;
}
inline void MetaData::clear_key() {
  if (key_ != &::google::protobuf::internal::kEmptyString) {
    key_->clear();
  }
  clear_has_key();
}
inline const ::std::string& MetaData::key() const {
  return *key_;
}
inline void MetaData::set_key(const ::std::string& value) {
  set_has_key();
  if (key_ == &::google::protobuf::internal::kEmptyString) {
    key_ = new ::std::string;
  }
  key_->assign(value);
}
inline void MetaData::set_key(const char* value) {
  set_has_key();
  if (key_ == &::google::protobuf::internal::kEmptyString) {
    key_ = new ::std::string;
  }
  key_->assign(value);
}
inline void MetaData::set_key(const char* value, size_t size) {
  set_has_key();
  if (key_ == &::google::protobuf::internal::kEmptyString) {
    key_ = new ::std::string;
  }
  key_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* MetaData::mutable_key() {
  set_has_key();
  if (key_ == &::google::protobuf::internal::kEmptyString) {
    key_ = new ::std::string;
  }
  return key_;
}
inline ::std::string* MetaData::release_key() {
  clear_has_key();
  if (key_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = key_;
    key_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}

// required string value = 2;
inline bool MetaData::has_value() const {
  return (_has_bits_[0] & 0x00000002u) != 0;
}
inline void MetaData::set_has_value() {
  _has_bits_[0] |= 0x00000002u;
}
inline void MetaData::clear_has_value() {
  _has_bits_[0] &= ~0x00000002u;
}
inline void MetaData::clear_value() {
  if (value_ != &::google::protobuf::internal::kEmptyString) {
    value_->clear();
  }
  clear_has_value();
}
inline const ::std::string& MetaData::value() const {
  return *value_;
}
inline void MetaData::set_value(const ::std::string& value) {
  set_has_value();
  if (value_ == &::google::protobuf::internal::kEmptyString) {
    value_ = new ::std::string;
  }
  value_->assign(value);
}
inline void MetaData::set_value(const char* value) {
  set_has_value();
  if (value_ == &::google::protobuf::internal::kEmptyString) {
    value_ = new ::std::string;
  }
  value_->assign(value);
}
inline void MetaData::set_value(const char* value, size_t size) {
  set_has_value();
  if (value_ == &::google::protobuf::internal::kEmptyString) {
    value_ = new ::std::string;
  }
  value_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* MetaData::mutable_value() {
  set_has_value();
  if (value_ == &::google::protobuf::internal::kEmptyString) {
    value_ = new ::std::string;
  }
  return value_;
}
inline ::std::string* MetaData::release_value() {
  clear_has_value();
  if (value_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = value_;
    value_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}


// @@protoc_insertion_point(namespace_scope)

}  // namespace goby

#ifndef SWIG
namespace google {
namespace protobuf {


}  // namespace google
}  // namespace protobuf
#endif  // SWIG

// @@protoc_insertion_point(global_scope)

#endif  // PROTOBUF_Reads_2eproto__INCLUDED
