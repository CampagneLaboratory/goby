#
# Environment variables used for SGE alignment jobs
#

#
# Email address to notify upon job completion
#
SGE_STATUS_MAILTO=

#
# Memory size for SGE
#
SGE_MEMORY=3g

#
# JVM flags (memory should match SGE_MEMORY)
#
JVM_FLAGS=-Xmx3g

#
# Aligner
#
ALIGNER=bwa

#
# Full path to native aligner executable
#
BWA_ALIGNER_EXEC=
LASTAG_ALIGNER_EXEC=

#
# Colorspace option (comment out or leave blank if not using colorspace)
#
#COLORSPACE=--color-space
COLORSPACE=

#
# Full path to reference file in compact format
#
REFERENCE=

#
# Full path to the index directory (i.e., reference cache)
#
REFERENCE_INDEX_DIRECTORY=

#
# Full path to the reads file
#
READS=

#
# Number of bytes to read per chunk
# If this is not set then the entire read is processed in a single job
#
CHUNK_SIZE=10000000
