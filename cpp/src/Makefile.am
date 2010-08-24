## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libgoby.la

AM_CXXFLAGS = @PROTOBUF_CFLAGS@ ${BOOST_CPPFLAGS}

libgoby_la_LIBADD = @PROTOBUF_LIBS@ ${BOOST_LDFLAGS}
libgoby_la_LDFLAGS = -version-info @version_info@ -export-dynamic # -no-undefined
libgoby_la_SOURCES = \
	goby/Alignments.cc \
	goby/Alignments.h \
	goby/Alignments.pb.cc \
	goby/Alignments.pb.h \
	goby/common.cc \
	goby/MessageChunks.h \
	goby/Reads.cc \
	goby/Reads.h \
	goby/Reads.pb.cc \
	goby/Reads.pb.h \
	goby/TooManyHits.cc \
	goby/TooManyHits.h

nobase_include_HEADERS = \
	goby/Alignments.h \
	goby/Alignments.pb.h \
	goby/common.h \
	goby/hash.h \
	goby/MessageChunks.h \
	goby/Reads.h \
	goby/Reads.pb.h \
	goby/TooManyHits.h

bin_PROGRAMS = GobyAlignmentStats GobyAlignmentToText GobyReadsStats GobyFastaToCompact

GobyAlignmentStats_LDADD = libgoby.la
GobyAlignmentStats_SOURCES = \
	GobyAlignmentStats.cc

GobyAlignmentToText_LDADD = libgoby.la
GobyAlignmentToText_SOURCES = \
	GobyAlignmentToText.cc

GobyReadsStats_LDADD = libgoby.la
GobyReadsStats_SOURCES = \
	GobyReadsStats.cc

GobyFastaToCompact_LDADD = libgoby.la ${BOOST_PROGRAM_OPTIONS_LIB}
GobyFastaToCompact_SOURCES = \
	GobyFastaToCompact.cc

MAINTAINERCLEANFILES = \
	Makefile.in
