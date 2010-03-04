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
SGE_MEMORY=4g

#
# JVM flags (memory should match SGE_MEMORY)
#
SGE_JVM_FLAGS=-Xmx4g

#
# Aligner
#
ALIGNER=bwa

#
# Full path to native aligner executable directory
#
BWA_ALIGNER_PATH=
LAST_ALIGNER_PATH=
LASTAG_ALIGNER_PATH=

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
# Full path to directory containing one or more transcripts in compact format
#
TRANSCRIPT_DIRECTORY=

#
# Full path to the index directory (i.e., reference cache) for transcripts
#
TRANSCRIPT_INDEX_DIRECTORY=

#
# Full path to gene to transcript mapping file
#
GENE_TRANSCRIPT_MAP_FILE=

#
# Number of bytes to read per chunk
# If this is not set then the entire read is processed in a single job
#
CHUNK_SIZE=10000000
