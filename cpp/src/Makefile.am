## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libgoby.la

BUILT_SOURCES = goby.timestamp

AM_CXXFLAGS = @PROTOBUF_CFLAGS@ @LIBPCRECPP_CFLAGS@ ${BOOST_CPPFLAGS}
AM_CFLAGS = @PROTOBUF_CFLAGS@ @LIBPCRECPP_CFLAGS@ ${BOOST_CPPFLAGS}

libgoby_la_LIBADD = @PROTOBUF_LIBS@ @LIBPCRECPP_LIBS@
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
	goby/TooManyHits.h \
	goby/C_Reads.cc \
	goby/C_Reads.h \
	goby/C_Alignments.cc \
	goby/C_Alignments.h \
	goby/C_CompactHelpers.h \
	goby/SamFlags.h

nobase_include_HEADERS = \
	goby/Alignments.h \
	goby/Alignments.pb.h \
	goby/common.h \
	goby/hash.h \
	goby/MessageChunks.h \
	goby/Reads.h \
	goby/Reads.pb.h \
	goby/TooManyHits.h \
	goby/C_Reads.h \
	goby/C_Alignments.h \
	goby/C_CompactHelpers.h \
	goby/SamFlags.h \
	goby/timestamp.h

bin_PROGRAMS = GobyAlignmentStats GobyAlignmentToText GobyReadsStats GobyFastaToCompact C_Sample_CompactReads

GobyAlignmentStats_LDADD = libgoby.la ${BOOST_LDFLAGS} ${BOOST_SYSTEM_LIB} ${BOOST_DATE_TIME_LIB} ${BOOST_FILESYSTEM_LIB}
GobyAlignmentStats_SOURCES = \
	GobyAlignmentStats.cc

GobyAlignmentToText_LDADD = libgoby.la ${BOOST_LDFLAGS} ${BOOST_SYSTEM_LIB} ${BOOST_DATE_TIME_LIB} ${BOOST_FILESYSTEM_LIB}
GobyAlignmentToText_SOURCES = \
	GobyAlignmentToText.cc

GobyReadsStats_LDADD = libgoby.la ${BOOST_LDFLAGS} ${BOOST_SYSTEM_LIB} ${BOOST_DATE_TIME_LIB} ${BOOST_FILESYSTEM_LIB}
GobyReadsStats_SOURCES = \
	GobyReadsStats.cc

GobyFastaToCompact_LDADD = libgoby.la ${BOOST_LDFLAGS} ${BOOST_SYSTEM_LIB} ${BOOST_DATE_TIME_LIB} ${BOOST_FILESYSTEM_LIB} ${BOOST_PROGRAM_OPTIONS_LIB}
GobyFastaToCompact_SOURCES = \
	GobyFastaToCompact.cc

C_Sample_CompactReads_LDADD = libgoby.la ${BOOST_LDFLAGS} ${BOOST_SYSTEM_LIB} ${BOOST_DATE_TIME_LIB} ${BOOST_FILESYSTEM_LIB} ${BOOST_PROGRAM_OPTIONS_LIB}
C_Sample_CompactReads_SOURCES = \
	C_Sample_CompactReads.c

#
# TODO:
# How to ALWAYS make this file at the top of the build?
# but nothing should "depend" on it, so just the change of this
# file should not trigger other things to build.
#
goby.timestamp:
	echo '#define BUILD_TIMESTAMP "'`date "+%G%m%d%k%M%S"`'"' > goby/timestamp.h

MAINTAINERCLEANFILES = \
	Makefile.in \
	goby/timestamp.h

CLEANFILES = \
	goby/timestamp.h


