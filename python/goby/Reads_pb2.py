# Generated by the protocol buffer compiler.  DO NOT EDIT!

from google.protobuf import descriptor
from google.protobuf import message
from google.protobuf import reflection
from google.protobuf import descriptor_pb2
# @@protoc_insertion_point(imports)


DESCRIPTOR = descriptor.FileDescriptor(
  name='Reads.proto',
  package='goby',
  serialized_pb='\n\x0bReads.proto\x12\x04goby\"0\n\x0eReadCollection\x12\x1e\n\x05reads\x18\x01 \x03(\x0b\x32\x0f.goby.ReadEntry\"\xf1\x01\n\tReadEntry\x12\x12\n\nread_index\x18\x01 \x02(\r\x12\x15\n\rbarcode_index\x18\n \x01(\r\x12\x17\n\x0fread_identifier\x18\x17 \x01(\t\x12\x13\n\x0b\x64\x65scription\x18\x16 \x01(\t\x12\x13\n\x0bread_length\x18\x02 \x02(\r\x12\x10\n\x08sequence\x18\x03 \x01(\x0c\x12\x15\n\rsequence_pair\x18\x05 \x01(\x0c\x12\x18\n\x10read_length_pair\x18\x06 \x01(\r\x12\x16\n\x0equality_scores\x18\x04 \x01(\x0c\x12\x1b\n\x13quality_scores_pair\x18\x07 \x01(\x0c\x42\"\n\x1e\x65\x64u.cornell.med.icb.goby.readsH\x01')




_READCOLLECTION = descriptor.Descriptor(
  name='ReadCollection',
  full_name='goby.ReadCollection',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    descriptor.FieldDescriptor(
      name='reads', full_name='goby.ReadCollection.reads', index=0,
      number=1, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  extension_ranges=[],
  serialized_start=21,
  serialized_end=69,
)


_READENTRY = descriptor.Descriptor(
  name='ReadEntry',
  full_name='goby.ReadEntry',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    descriptor.FieldDescriptor(
      name='read_index', full_name='goby.ReadEntry.read_index', index=0,
      number=1, type=13, cpp_type=3, label=2,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    descriptor.FieldDescriptor(
      name='barcode_index', full_name='goby.ReadEntry.barcode_index', index=1,
      number=10, type=13, cpp_type=3, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    descriptor.FieldDescriptor(
      name='read_identifier', full_name='goby.ReadEntry.read_identifier', index=2,
      number=23, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=unicode("", "utf-8"),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    descriptor.FieldDescriptor(
      name='description', full_name='goby.ReadEntry.description', index=3,
      number=22, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=unicode("", "utf-8"),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    descriptor.FieldDescriptor(
      name='read_length', full_name='goby.ReadEntry.read_length', index=4,
      number=2, type=13, cpp_type=3, label=2,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    descriptor.FieldDescriptor(
      name='sequence', full_name='goby.ReadEntry.sequence', index=5,
      number=3, type=12, cpp_type=9, label=1,
      has_default_value=False, default_value="",
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    descriptor.FieldDescriptor(
      name='sequence_pair', full_name='goby.ReadEntry.sequence_pair', index=6,
      number=5, type=12, cpp_type=9, label=1,
      has_default_value=False, default_value="",
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    descriptor.FieldDescriptor(
      name='read_length_pair', full_name='goby.ReadEntry.read_length_pair', index=7,
      number=6, type=13, cpp_type=3, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    descriptor.FieldDescriptor(
      name='quality_scores', full_name='goby.ReadEntry.quality_scores', index=8,
      number=4, type=12, cpp_type=9, label=1,
      has_default_value=False, default_value="",
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    descriptor.FieldDescriptor(
      name='quality_scores_pair', full_name='goby.ReadEntry.quality_scores_pair', index=9,
      number=7, type=12, cpp_type=9, label=1,
      has_default_value=False, default_value="",
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  extension_ranges=[],
  serialized_start=72,
  serialized_end=313,
)


_READCOLLECTION.fields_by_name['reads'].message_type = _READENTRY

class ReadCollection(message.Message):
  __metaclass__ = reflection.GeneratedProtocolMessageType
  DESCRIPTOR = _READCOLLECTION
  
  # @@protoc_insertion_point(class_scope:goby.ReadCollection)

class ReadEntry(message.Message):
  __metaclass__ = reflection.GeneratedProtocolMessageType
  DESCRIPTOR = _READENTRY
  
  # @@protoc_insertion_point(class_scope:goby.ReadEntry)

# @@protoc_insertion_point(module_scope)
