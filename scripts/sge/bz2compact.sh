#!/bin/sh

#
# Script to convert a fastq.bz2 file to compact format. The second argument must provide the quality encoding
# used in the fastq file.
# Example of use:
# foreach file (full-path/*.fastq.bz2)
#  bz2compact.sh $file (SANGER|ILLUMINA|SOLEXA)
# end
# The previous lines will convert each file with extension .fastq.bz2 to the compact-reads format, running
# in parallel on an SGE grid. The file bz2compact.sh must be on the path.


export WORKING_DIR=`dirname "$0"`
export GOBY=`which goby`
export BASENAME=`basename $1 .fastq.bz2`
export DATA_DIR=`dirname $1`
export QUALITY_ENCODING=$2

CMD=`qsub -terse -v GOBY=$GOBY -v BASENAME=$BASENAME -v DATA_DIR=$DATA_DIR -v QUALITY_ENCODING=$QUALITY_ENCODING ${WORKING_DIR}/bz2compact.qsub `
echo $CMD


