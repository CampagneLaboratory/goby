// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: src/edu/cornell/med/icb/goby/reads/Reads.proto

package edu.cornell.med.icb.goby.reads;

public final class Reads {
  private Reads() {}
  public static void registerAllExtensions(
      com.google.protobuf.ExtensionRegistry registry) {
  }
  public static final class ReadCollection extends
      com.google.protobuf.GeneratedMessage {
    // Use ReadCollection.newBuilder() to construct.
    private ReadCollection() {
      initFields();
    }
    private ReadCollection(boolean noInit) {}
    
    private static final ReadCollection defaultInstance;
    public static ReadCollection getDefaultInstance() {
      return defaultInstance;
    }
    
    public ReadCollection getDefaultInstanceForType() {
      return defaultInstance;
    }
    
    public static final com.google.protobuf.Descriptors.Descriptor
        getDescriptor() {
      return edu.cornell.med.icb.goby.reads.Reads.internal_static_goby_ReadCollection_descriptor;
    }
    
    protected com.google.protobuf.GeneratedMessage.FieldAccessorTable
        internalGetFieldAccessorTable() {
      return edu.cornell.med.icb.goby.reads.Reads.internal_static_goby_ReadCollection_fieldAccessorTable;
    }
    
    // repeated .goby.ReadEntry reads = 1;
    public static final int READS_FIELD_NUMBER = 1;
    private java.util.List<edu.cornell.med.icb.goby.reads.Reads.ReadEntry> reads_ =
      java.util.Collections.emptyList();
    public java.util.List<edu.cornell.med.icb.goby.reads.Reads.ReadEntry> getReadsList() {
      return reads_;
    }
    public int getReadsCount() { return reads_.size(); }
    public edu.cornell.med.icb.goby.reads.Reads.ReadEntry getReads(int index) {
      return reads_.get(index);
    }
    
    private void initFields() {
    }
    public final boolean isInitialized() {
      for (edu.cornell.med.icb.goby.reads.Reads.ReadEntry element : getReadsList()) {
        if (!element.isInitialized()) return false;
      }
      return true;
    }
    
    public void writeTo(com.google.protobuf.CodedOutputStream output)
                        throws java.io.IOException {
      getSerializedSize();
      for (edu.cornell.med.icb.goby.reads.Reads.ReadEntry element : getReadsList()) {
        output.writeMessage(1, element);
      }
      getUnknownFields().writeTo(output);
    }
    
    private int memoizedSerializedSize = -1;
    public int getSerializedSize() {
      int size = memoizedSerializedSize;
      if (size != -1) return size;
    
      size = 0;
      for (edu.cornell.med.icb.goby.reads.Reads.ReadEntry element : getReadsList()) {
        size += com.google.protobuf.CodedOutputStream
          .computeMessageSize(1, element);
      }
      size += getUnknownFields().getSerializedSize();
      memoizedSerializedSize = size;
      return size;
    }
    
    public static edu.cornell.med.icb.goby.reads.Reads.ReadCollection parseFrom(
        com.google.protobuf.ByteString data)
        throws com.google.protobuf.InvalidProtocolBufferException {
      return newBuilder().mergeFrom(data).buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadCollection parseFrom(
        com.google.protobuf.ByteString data,
        com.google.protobuf.ExtensionRegistryLite extensionRegistry)
        throws com.google.protobuf.InvalidProtocolBufferException {
      return newBuilder().mergeFrom(data, extensionRegistry)
               .buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadCollection parseFrom(byte[] data)
        throws com.google.protobuf.InvalidProtocolBufferException {
      return newBuilder().mergeFrom(data).buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadCollection parseFrom(
        byte[] data,
        com.google.protobuf.ExtensionRegistryLite extensionRegistry)
        throws com.google.protobuf.InvalidProtocolBufferException {
      return newBuilder().mergeFrom(data, extensionRegistry)
               .buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadCollection parseFrom(java.io.InputStream input)
        throws java.io.IOException {
      return newBuilder().mergeFrom(input).buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadCollection parseFrom(
        java.io.InputStream input,
        com.google.protobuf.ExtensionRegistryLite extensionRegistry)
        throws java.io.IOException {
      return newBuilder().mergeFrom(input, extensionRegistry)
               .buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadCollection parseDelimitedFrom(java.io.InputStream input)
        throws java.io.IOException {
      Builder builder = newBuilder();
      if (builder.mergeDelimitedFrom(input)) {
        return builder.buildParsed();
      } else {
        return null;
      }
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadCollection parseDelimitedFrom(
        java.io.InputStream input,
        com.google.protobuf.ExtensionRegistryLite extensionRegistry)
        throws java.io.IOException {
      Builder builder = newBuilder();
      if (builder.mergeDelimitedFrom(input, extensionRegistry)) {
        return builder.buildParsed();
      } else {
        return null;
      }
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadCollection parseFrom(
        com.google.protobuf.CodedInputStream input)
        throws java.io.IOException {
      return newBuilder().mergeFrom(input).buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadCollection parseFrom(
        com.google.protobuf.CodedInputStream input,
        com.google.protobuf.ExtensionRegistryLite extensionRegistry)
        throws java.io.IOException {
      return newBuilder().mergeFrom(input, extensionRegistry)
               .buildParsed();
    }
    
    public static Builder newBuilder() { return Builder.create(); }
    public Builder newBuilderForType() { return newBuilder(); }
    public static Builder newBuilder(edu.cornell.med.icb.goby.reads.Reads.ReadCollection prototype) {
      return newBuilder().mergeFrom(prototype);
    }
    public Builder toBuilder() { return newBuilder(this); }
    
    public static final class Builder extends
        com.google.protobuf.GeneratedMessage.Builder<Builder> {
      private edu.cornell.med.icb.goby.reads.Reads.ReadCollection result;
      
      // Construct using edu.cornell.med.icb.goby.reads.Reads.ReadCollection.newBuilder()
      private Builder() {}
      
      private static Builder create() {
        Builder builder = new Builder();
        builder.result = new edu.cornell.med.icb.goby.reads.Reads.ReadCollection();
        return builder;
      }
      
      protected edu.cornell.med.icb.goby.reads.Reads.ReadCollection internalGetResult() {
        return result;
      }
      
      public Builder clear() {
        if (result == null) {
          throw new IllegalStateException(
            "Cannot call clear() after build().");
        }
        result = new edu.cornell.med.icb.goby.reads.Reads.ReadCollection();
        return this;
      }
      
      public Builder clone() {
        return create().mergeFrom(result);
      }
      
      public com.google.protobuf.Descriptors.Descriptor
          getDescriptorForType() {
        return edu.cornell.med.icb.goby.reads.Reads.ReadCollection.getDescriptor();
      }
      
      public edu.cornell.med.icb.goby.reads.Reads.ReadCollection getDefaultInstanceForType() {
        return edu.cornell.med.icb.goby.reads.Reads.ReadCollection.getDefaultInstance();
      }
      
      public boolean isInitialized() {
        return result.isInitialized();
      }
      public edu.cornell.med.icb.goby.reads.Reads.ReadCollection build() {
        if (result != null && !isInitialized()) {
          throw newUninitializedMessageException(result);
        }
        return buildPartial();
      }
      
      private edu.cornell.med.icb.goby.reads.Reads.ReadCollection buildParsed()
          throws com.google.protobuf.InvalidProtocolBufferException {
        if (!isInitialized()) {
          throw newUninitializedMessageException(
            result).asInvalidProtocolBufferException();
        }
        return buildPartial();
      }
      
      public edu.cornell.med.icb.goby.reads.Reads.ReadCollection buildPartial() {
        if (result == null) {
          throw new IllegalStateException(
            "build() has already been called on this Builder.");
        }
        if (result.reads_ != java.util.Collections.EMPTY_LIST) {
          result.reads_ =
            java.util.Collections.unmodifiableList(result.reads_);
        }
        edu.cornell.med.icb.goby.reads.Reads.ReadCollection returnMe = result;
        result = null;
        return returnMe;
      }
      
      public Builder mergeFrom(com.google.protobuf.Message other) {
        if (other instanceof edu.cornell.med.icb.goby.reads.Reads.ReadCollection) {
          return mergeFrom((edu.cornell.med.icb.goby.reads.Reads.ReadCollection)other);
        } else {
          super.mergeFrom(other);
          return this;
        }
      }
      
      public Builder mergeFrom(edu.cornell.med.icb.goby.reads.Reads.ReadCollection other) {
        if (other == edu.cornell.med.icb.goby.reads.Reads.ReadCollection.getDefaultInstance()) return this;
        if (!other.reads_.isEmpty()) {
          if (result.reads_.isEmpty()) {
            result.reads_ = new java.util.ArrayList<edu.cornell.med.icb.goby.reads.Reads.ReadEntry>();
          }
          result.reads_.addAll(other.reads_);
        }
        this.mergeUnknownFields(other.getUnknownFields());
        return this;
      }
      
      public Builder mergeFrom(
          com.google.protobuf.CodedInputStream input,
          com.google.protobuf.ExtensionRegistryLite extensionRegistry)
          throws java.io.IOException {
        com.google.protobuf.UnknownFieldSet.Builder unknownFields =
          com.google.protobuf.UnknownFieldSet.newBuilder(
            this.getUnknownFields());
        while (true) {
          int tag = input.readTag();
          switch (tag) {
            case 0:
              this.setUnknownFields(unknownFields.build());
              return this;
            default: {
              if (!parseUnknownField(input, unknownFields,
                                     extensionRegistry, tag)) {
                this.setUnknownFields(unknownFields.build());
                return this;
              }
              break;
            }
            case 10: {
              edu.cornell.med.icb.goby.reads.Reads.ReadEntry.Builder subBuilder = edu.cornell.med.icb.goby.reads.Reads.ReadEntry.newBuilder();
              input.readMessage(subBuilder, extensionRegistry);
              addReads(subBuilder.buildPartial());
              break;
            }
          }
        }
      }
      
      
      // repeated .goby.ReadEntry reads = 1;
      public java.util.List<edu.cornell.med.icb.goby.reads.Reads.ReadEntry> getReadsList() {
        return java.util.Collections.unmodifiableList(result.reads_);
      }
      public int getReadsCount() {
        return result.getReadsCount();
      }
      public edu.cornell.med.icb.goby.reads.Reads.ReadEntry getReads(int index) {
        return result.getReads(index);
      }
      public Builder setReads(int index, edu.cornell.med.icb.goby.reads.Reads.ReadEntry value) {
        if (value == null) {
          throw new NullPointerException();
        }
        result.reads_.set(index, value);
        return this;
      }
      public Builder setReads(int index, edu.cornell.med.icb.goby.reads.Reads.ReadEntry.Builder builderForValue) {
        result.reads_.set(index, builderForValue.build());
        return this;
      }
      public Builder addReads(edu.cornell.med.icb.goby.reads.Reads.ReadEntry value) {
        if (value == null) {
          throw new NullPointerException();
        }
        if (result.reads_.isEmpty()) {
          result.reads_ = new java.util.ArrayList<edu.cornell.med.icb.goby.reads.Reads.ReadEntry>();
        }
        result.reads_.add(value);
        return this;
      }
      public Builder addReads(edu.cornell.med.icb.goby.reads.Reads.ReadEntry.Builder builderForValue) {
        if (result.reads_.isEmpty()) {
          result.reads_ = new java.util.ArrayList<edu.cornell.med.icb.goby.reads.Reads.ReadEntry>();
        }
        result.reads_.add(builderForValue.build());
        return this;
      }
      public Builder addAllReads(
          java.lang.Iterable<? extends edu.cornell.med.icb.goby.reads.Reads.ReadEntry> values) {
        if (result.reads_.isEmpty()) {
          result.reads_ = new java.util.ArrayList<edu.cornell.med.icb.goby.reads.Reads.ReadEntry>();
        }
        super.addAll(values, result.reads_);
        return this;
      }
      public Builder clearReads() {
        result.reads_ = java.util.Collections.emptyList();
        return this;
      }
      
      // @@protoc_insertion_point(builder_scope:goby.ReadCollection)
    }
    
    static {
      defaultInstance = new ReadCollection(true);
      edu.cornell.med.icb.goby.reads.Reads.internalForceInit();
      defaultInstance.initFields();
    }
    
    // @@protoc_insertion_point(class_scope:goby.ReadCollection)
  }
  
  public static final class ReadEntry extends
      com.google.protobuf.GeneratedMessage {
    // Use ReadEntry.newBuilder() to construct.
    private ReadEntry() {
      initFields();
    }
    private ReadEntry(boolean noInit) {}
    
    private static final ReadEntry defaultInstance;
    public static ReadEntry getDefaultInstance() {
      return defaultInstance;
    }
    
    public ReadEntry getDefaultInstanceForType() {
      return defaultInstance;
    }
    
    public static final com.google.protobuf.Descriptors.Descriptor
        getDescriptor() {
      return edu.cornell.med.icb.goby.reads.Reads.internal_static_goby_ReadEntry_descriptor;
    }
    
    protected com.google.protobuf.GeneratedMessage.FieldAccessorTable
        internalGetFieldAccessorTable() {
      return edu.cornell.med.icb.goby.reads.Reads.internal_static_goby_ReadEntry_fieldAccessorTable;
    }
    
    // required uint32 readIndex = 1;
    public static final int READINDEX_FIELD_NUMBER = 1;
    private boolean hasReadIndex;
    private int readIndex_ = 0;
    public boolean hasReadIndex() { return hasReadIndex; }
    public int getReadIndex() { return readIndex_; }
    
    // optional uint32 barcodeIndex = 10;
    public static final int BARCODEINDEX_FIELD_NUMBER = 10;
    private boolean hasBarcodeIndex;
    private int barcodeIndex_ = 0;
    public boolean hasBarcodeIndex() { return hasBarcodeIndex; }
    public int getBarcodeIndex() { return barcodeIndex_; }
    
    // optional string readIdentifier = 23;
    public static final int READIDENTIFIER_FIELD_NUMBER = 23;
    private boolean hasReadIdentifier;
    private java.lang.String readIdentifier_ = "";
    public boolean hasReadIdentifier() { return hasReadIdentifier; }
    public java.lang.String getReadIdentifier() { return readIdentifier_; }
    
    // optional string description = 22;
    public static final int DESCRIPTION_FIELD_NUMBER = 22;
    private boolean hasDescription;
    private java.lang.String description_ = "";
    public boolean hasDescription() { return hasDescription; }
    public java.lang.String getDescription() { return description_; }
    
    // required uint32 readLength = 2;
    public static final int READLENGTH_FIELD_NUMBER = 2;
    private boolean hasReadLength;
    private int readLength_ = 0;
    public boolean hasReadLength() { return hasReadLength; }
    public int getReadLength() { return readLength_; }
    
    // optional bytes sequence = 3;
    public static final int SEQUENCE_FIELD_NUMBER = 3;
    private boolean hasSequence;
    private com.google.protobuf.ByteString sequence_ = com.google.protobuf.ByteString.EMPTY;
    public boolean hasSequence() { return hasSequence; }
    public com.google.protobuf.ByteString getSequence() { return sequence_; }
    
    // optional bytes qualityScores = 4;
    public static final int QUALITYSCORES_FIELD_NUMBER = 4;
    private boolean hasQualityScores;
    private com.google.protobuf.ByteString qualityScores_ = com.google.protobuf.ByteString.EMPTY;
    public boolean hasQualityScores() { return hasQualityScores; }
    public com.google.protobuf.ByteString getQualityScores() { return qualityScores_; }
    
    private void initFields() {
    }
    public final boolean isInitialized() {
      if (!hasReadIndex) return false;
      if (!hasReadLength) return false;
      return true;
    }
    
    public void writeTo(com.google.protobuf.CodedOutputStream output)
                        throws java.io.IOException {
      getSerializedSize();
      if (hasReadIndex()) {
        output.writeUInt32(1, getReadIndex());
      }
      if (hasReadLength()) {
        output.writeUInt32(2, getReadLength());
      }
      if (hasSequence()) {
        output.writeBytes(3, getSequence());
      }
      if (hasQualityScores()) {
        output.writeBytes(4, getQualityScores());
      }
      if (hasBarcodeIndex()) {
        output.writeUInt32(10, getBarcodeIndex());
      }
      if (hasDescription()) {
        output.writeString(22, getDescription());
      }
      if (hasReadIdentifier()) {
        output.writeString(23, getReadIdentifier());
      }
      getUnknownFields().writeTo(output);
    }
    
    private int memoizedSerializedSize = -1;
    public int getSerializedSize() {
      int size = memoizedSerializedSize;
      if (size != -1) return size;
    
      size = 0;
      if (hasReadIndex()) {
        size += com.google.protobuf.CodedOutputStream
          .computeUInt32Size(1, getReadIndex());
      }
      if (hasReadLength()) {
        size += com.google.protobuf.CodedOutputStream
          .computeUInt32Size(2, getReadLength());
      }
      if (hasSequence()) {
        size += com.google.protobuf.CodedOutputStream
          .computeBytesSize(3, getSequence());
      }
      if (hasQualityScores()) {
        size += com.google.protobuf.CodedOutputStream
          .computeBytesSize(4, getQualityScores());
      }
      if (hasBarcodeIndex()) {
        size += com.google.protobuf.CodedOutputStream
          .computeUInt32Size(10, getBarcodeIndex());
      }
      if (hasDescription()) {
        size += com.google.protobuf.CodedOutputStream
          .computeStringSize(22, getDescription());
      }
      if (hasReadIdentifier()) {
        size += com.google.protobuf.CodedOutputStream
          .computeStringSize(23, getReadIdentifier());
      }
      size += getUnknownFields().getSerializedSize();
      memoizedSerializedSize = size;
      return size;
    }
    
    public static edu.cornell.med.icb.goby.reads.Reads.ReadEntry parseFrom(
        com.google.protobuf.ByteString data)
        throws com.google.protobuf.InvalidProtocolBufferException {
      return newBuilder().mergeFrom(data).buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadEntry parseFrom(
        com.google.protobuf.ByteString data,
        com.google.protobuf.ExtensionRegistryLite extensionRegistry)
        throws com.google.protobuf.InvalidProtocolBufferException {
      return newBuilder().mergeFrom(data, extensionRegistry)
               .buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadEntry parseFrom(byte[] data)
        throws com.google.protobuf.InvalidProtocolBufferException {
      return newBuilder().mergeFrom(data).buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadEntry parseFrom(
        byte[] data,
        com.google.protobuf.ExtensionRegistryLite extensionRegistry)
        throws com.google.protobuf.InvalidProtocolBufferException {
      return newBuilder().mergeFrom(data, extensionRegistry)
               .buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadEntry parseFrom(java.io.InputStream input)
        throws java.io.IOException {
      return newBuilder().mergeFrom(input).buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadEntry parseFrom(
        java.io.InputStream input,
        com.google.protobuf.ExtensionRegistryLite extensionRegistry)
        throws java.io.IOException {
      return newBuilder().mergeFrom(input, extensionRegistry)
               .buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadEntry parseDelimitedFrom(java.io.InputStream input)
        throws java.io.IOException {
      Builder builder = newBuilder();
      if (builder.mergeDelimitedFrom(input)) {
        return builder.buildParsed();
      } else {
        return null;
      }
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadEntry parseDelimitedFrom(
        java.io.InputStream input,
        com.google.protobuf.ExtensionRegistryLite extensionRegistry)
        throws java.io.IOException {
      Builder builder = newBuilder();
      if (builder.mergeDelimitedFrom(input, extensionRegistry)) {
        return builder.buildParsed();
      } else {
        return null;
      }
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadEntry parseFrom(
        com.google.protobuf.CodedInputStream input)
        throws java.io.IOException {
      return newBuilder().mergeFrom(input).buildParsed();
    }
    public static edu.cornell.med.icb.goby.reads.Reads.ReadEntry parseFrom(
        com.google.protobuf.CodedInputStream input,
        com.google.protobuf.ExtensionRegistryLite extensionRegistry)
        throws java.io.IOException {
      return newBuilder().mergeFrom(input, extensionRegistry)
               .buildParsed();
    }
    
    public static Builder newBuilder() { return Builder.create(); }
    public Builder newBuilderForType() { return newBuilder(); }
    public static Builder newBuilder(edu.cornell.med.icb.goby.reads.Reads.ReadEntry prototype) {
      return newBuilder().mergeFrom(prototype);
    }
    public Builder toBuilder() { return newBuilder(this); }
    
    public static final class Builder extends
        com.google.protobuf.GeneratedMessage.Builder<Builder> {
      private edu.cornell.med.icb.goby.reads.Reads.ReadEntry result;
      
      // Construct using edu.cornell.med.icb.goby.reads.Reads.ReadEntry.newBuilder()
      private Builder() {}
      
      private static Builder create() {
        Builder builder = new Builder();
        builder.result = new edu.cornell.med.icb.goby.reads.Reads.ReadEntry();
        return builder;
      }
      
      protected edu.cornell.med.icb.goby.reads.Reads.ReadEntry internalGetResult() {
        return result;
      }
      
      public Builder clear() {
        if (result == null) {
          throw new IllegalStateException(
            "Cannot call clear() after build().");
        }
        result = new edu.cornell.med.icb.goby.reads.Reads.ReadEntry();
        return this;
      }
      
      public Builder clone() {
        return create().mergeFrom(result);
      }
      
      public com.google.protobuf.Descriptors.Descriptor
          getDescriptorForType() {
        return edu.cornell.med.icb.goby.reads.Reads.ReadEntry.getDescriptor();
      }
      
      public edu.cornell.med.icb.goby.reads.Reads.ReadEntry getDefaultInstanceForType() {
        return edu.cornell.med.icb.goby.reads.Reads.ReadEntry.getDefaultInstance();
      }
      
      public boolean isInitialized() {
        return result.isInitialized();
      }
      public edu.cornell.med.icb.goby.reads.Reads.ReadEntry build() {
        if (result != null && !isInitialized()) {
          throw newUninitializedMessageException(result);
        }
        return buildPartial();
      }
      
      private edu.cornell.med.icb.goby.reads.Reads.ReadEntry buildParsed()
          throws com.google.protobuf.InvalidProtocolBufferException {
        if (!isInitialized()) {
          throw newUninitializedMessageException(
            result).asInvalidProtocolBufferException();
        }
        return buildPartial();
      }
      
      public edu.cornell.med.icb.goby.reads.Reads.ReadEntry buildPartial() {
        if (result == null) {
          throw new IllegalStateException(
            "build() has already been called on this Builder.");
        }
        edu.cornell.med.icb.goby.reads.Reads.ReadEntry returnMe = result;
        result = null;
        return returnMe;
      }
      
      public Builder mergeFrom(com.google.protobuf.Message other) {
        if (other instanceof edu.cornell.med.icb.goby.reads.Reads.ReadEntry) {
          return mergeFrom((edu.cornell.med.icb.goby.reads.Reads.ReadEntry)other);
        } else {
          super.mergeFrom(other);
          return this;
        }
      }
      
      public Builder mergeFrom(edu.cornell.med.icb.goby.reads.Reads.ReadEntry other) {
        if (other == edu.cornell.med.icb.goby.reads.Reads.ReadEntry.getDefaultInstance()) return this;
        if (other.hasReadIndex()) {
          setReadIndex(other.getReadIndex());
        }
        if (other.hasBarcodeIndex()) {
          setBarcodeIndex(other.getBarcodeIndex());
        }
        if (other.hasReadIdentifier()) {
          setReadIdentifier(other.getReadIdentifier());
        }
        if (other.hasDescription()) {
          setDescription(other.getDescription());
        }
        if (other.hasReadLength()) {
          setReadLength(other.getReadLength());
        }
        if (other.hasSequence()) {
          setSequence(other.getSequence());
        }
        if (other.hasQualityScores()) {
          setQualityScores(other.getQualityScores());
        }
        this.mergeUnknownFields(other.getUnknownFields());
        return this;
      }
      
      public Builder mergeFrom(
          com.google.protobuf.CodedInputStream input,
          com.google.protobuf.ExtensionRegistryLite extensionRegistry)
          throws java.io.IOException {
        com.google.protobuf.UnknownFieldSet.Builder unknownFields =
          com.google.protobuf.UnknownFieldSet.newBuilder(
            this.getUnknownFields());
        while (true) {
          int tag = input.readTag();
          switch (tag) {
            case 0:
              this.setUnknownFields(unknownFields.build());
              return this;
            default: {
              if (!parseUnknownField(input, unknownFields,
                                     extensionRegistry, tag)) {
                this.setUnknownFields(unknownFields.build());
                return this;
              }
              break;
            }
            case 8: {
              setReadIndex(input.readUInt32());
              break;
            }
            case 16: {
              setReadLength(input.readUInt32());
              break;
            }
            case 26: {
              setSequence(input.readBytes());
              break;
            }
            case 34: {
              setQualityScores(input.readBytes());
              break;
            }
            case 80: {
              setBarcodeIndex(input.readUInt32());
              break;
            }
            case 178: {
              setDescription(input.readString());
              break;
            }
            case 186: {
              setReadIdentifier(input.readString());
              break;
            }
          }
        }
      }
      
      
      // required uint32 readIndex = 1;
      public boolean hasReadIndex() {
        return result.hasReadIndex();
      }
      public int getReadIndex() {
        return result.getReadIndex();
      }
      public Builder setReadIndex(int value) {
        result.hasReadIndex = true;
        result.readIndex_ = value;
        return this;
      }
      public Builder clearReadIndex() {
        result.hasReadIndex = false;
        result.readIndex_ = 0;
        return this;
      }
      
      // optional uint32 barcodeIndex = 10;
      public boolean hasBarcodeIndex() {
        return result.hasBarcodeIndex();
      }
      public int getBarcodeIndex() {
        return result.getBarcodeIndex();
      }
      public Builder setBarcodeIndex(int value) {
        result.hasBarcodeIndex = true;
        result.barcodeIndex_ = value;
        return this;
      }
      public Builder clearBarcodeIndex() {
        result.hasBarcodeIndex = false;
        result.barcodeIndex_ = 0;
        return this;
      }
      
      // optional string readIdentifier = 23;
      public boolean hasReadIdentifier() {
        return result.hasReadIdentifier();
      }
      public java.lang.String getReadIdentifier() {
        return result.getReadIdentifier();
      }
      public Builder setReadIdentifier(java.lang.String value) {
        if (value == null) {
    throw new NullPointerException();
  }
  result.hasReadIdentifier = true;
        result.readIdentifier_ = value;
        return this;
      }
      public Builder clearReadIdentifier() {
        result.hasReadIdentifier = false;
        result.readIdentifier_ = getDefaultInstance().getReadIdentifier();
        return this;
      }
      
      // optional string description = 22;
      public boolean hasDescription() {
        return result.hasDescription();
      }
      public java.lang.String getDescription() {
        return result.getDescription();
      }
      public Builder setDescription(java.lang.String value) {
        if (value == null) {
    throw new NullPointerException();
  }
  result.hasDescription = true;
        result.description_ = value;
        return this;
      }
      public Builder clearDescription() {
        result.hasDescription = false;
        result.description_ = getDefaultInstance().getDescription();
        return this;
      }
      
      // required uint32 readLength = 2;
      public boolean hasReadLength() {
        return result.hasReadLength();
      }
      public int getReadLength() {
        return result.getReadLength();
      }
      public Builder setReadLength(int value) {
        result.hasReadLength = true;
        result.readLength_ = value;
        return this;
      }
      public Builder clearReadLength() {
        result.hasReadLength = false;
        result.readLength_ = 0;
        return this;
      }
      
      // optional bytes sequence = 3;
      public boolean hasSequence() {
        return result.hasSequence();
      }
      public com.google.protobuf.ByteString getSequence() {
        return result.getSequence();
      }
      public Builder setSequence(com.google.protobuf.ByteString value) {
        if (value == null) {
    throw new NullPointerException();
  }
  result.hasSequence = true;
        result.sequence_ = value;
        return this;
      }
      public Builder clearSequence() {
        result.hasSequence = false;
        result.sequence_ = getDefaultInstance().getSequence();
        return this;
      }
      
      // optional bytes qualityScores = 4;
      public boolean hasQualityScores() {
        return result.hasQualityScores();
      }
      public com.google.protobuf.ByteString getQualityScores() {
        return result.getQualityScores();
      }
      public Builder setQualityScores(com.google.protobuf.ByteString value) {
        if (value == null) {
    throw new NullPointerException();
  }
  result.hasQualityScores = true;
        result.qualityScores_ = value;
        return this;
      }
      public Builder clearQualityScores() {
        result.hasQualityScores = false;
        result.qualityScores_ = getDefaultInstance().getQualityScores();
        return this;
      }
      
      // @@protoc_insertion_point(builder_scope:goby.ReadEntry)
    }
    
    static {
      defaultInstance = new ReadEntry(true);
      edu.cornell.med.icb.goby.reads.Reads.internalForceInit();
      defaultInstance.initFields();
    }
    
    // @@protoc_insertion_point(class_scope:goby.ReadEntry)
  }
  
  private static com.google.protobuf.Descriptors.Descriptor
    internal_static_goby_ReadCollection_descriptor;
  private static
    com.google.protobuf.GeneratedMessage.FieldAccessorTable
      internal_static_goby_ReadCollection_fieldAccessorTable;
  private static com.google.protobuf.Descriptors.Descriptor
    internal_static_goby_ReadEntry_descriptor;
  private static
    com.google.protobuf.GeneratedMessage.FieldAccessorTable
      internal_static_goby_ReadEntry_fieldAccessorTable;
  
  public static com.google.protobuf.Descriptors.FileDescriptor
      getDescriptor() {
    return descriptor;
  }
  private static com.google.protobuf.Descriptors.FileDescriptor
      descriptor;
  static {
    java.lang.String[] descriptorData = {
      "\n.src/edu/cornell/med/icb/goby/reads/Rea" +
      "ds.proto\022\004goby\"0\n\016ReadCollection\022\036\n\005read" +
      "s\030\001 \003(\0132\017.goby.ReadEntry\"\236\001\n\tReadEntry\022\021" +
      "\n\treadIndex\030\001 \002(\r\022\024\n\014barcodeIndex\030\n \001(\r\022" +
      "\026\n\016readIdentifier\030\027 \001(\t\022\023\n\013description\030\026" +
      " \001(\t\022\022\n\nreadLength\030\002 \002(\r\022\020\n\010sequence\030\003 \001" +
      "(\014\022\025\n\rqualityScores\030\004 \001(\014B\"\n\036edu.cornell" +
      ".med.icb.goby.readsH\001"
    };
    com.google.protobuf.Descriptors.FileDescriptor.InternalDescriptorAssigner assigner =
      new com.google.protobuf.Descriptors.FileDescriptor.InternalDescriptorAssigner() {
        public com.google.protobuf.ExtensionRegistry assignDescriptors(
            com.google.protobuf.Descriptors.FileDescriptor root) {
          descriptor = root;
          internal_static_goby_ReadCollection_descriptor =
            getDescriptor().getMessageTypes().get(0);
          internal_static_goby_ReadCollection_fieldAccessorTable = new
            com.google.protobuf.GeneratedMessage.FieldAccessorTable(
              internal_static_goby_ReadCollection_descriptor,
              new java.lang.String[] { "Reads", },
              edu.cornell.med.icb.goby.reads.Reads.ReadCollection.class,
              edu.cornell.med.icb.goby.reads.Reads.ReadCollection.Builder.class);
          internal_static_goby_ReadEntry_descriptor =
            getDescriptor().getMessageTypes().get(1);
          internal_static_goby_ReadEntry_fieldAccessorTable = new
            com.google.protobuf.GeneratedMessage.FieldAccessorTable(
              internal_static_goby_ReadEntry_descriptor,
              new java.lang.String[] { "ReadIndex", "BarcodeIndex", "ReadIdentifier", "Description", "ReadLength", "Sequence", "QualityScores", },
              edu.cornell.med.icb.goby.reads.Reads.ReadEntry.class,
              edu.cornell.med.icb.goby.reads.Reads.ReadEntry.Builder.class);
          return null;
        }
      };
    com.google.protobuf.Descriptors.FileDescriptor
      .internalBuildGeneratedFileFrom(descriptorData,
        new com.google.protobuf.Descriptors.FileDescriptor[] {
        }, assigner);
  }
  
  public static void internalForceInit() {}
  
  // @@protoc_insertion_point(outer_class_scope)
}
