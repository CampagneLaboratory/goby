## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libgoby.la

libgoby_la_CPPFLAGS = @PROTOBUF_CFLAGS@
libgoby_la_LIBADD = @PROTOBUF_LIBS@
libgoby_la_LDFLAGS = -version-info @version_info@ -export-dynamic -no-undefined
libgoby_la_SOURCES = \
	goby/Alignments.pb.cc \
	goby/Alignments.pb.h \
	goby/Reads.pb.cc \
	goby/Reads.pb.h

#bin_PROGRAMS = goby
#goby_SOURCES = \
	Alignments.pb.cc \
	Reads.pb.cc

MAINTAINERCLEANFILES = \
	Makefile.in
