#!/bin/sh

#
# Script to bundle goby and required data into a
# suitable format for submission to a PBS queue
#

# Absolute path to this script.
SCRIPT=$(readlink -f $0)
# Absolute path this script is in.
SCRIPT_DIR=$(dirname $SCRIPT)

if [ -z $1 ]; then
    echo "Job name is required"
    exit 1
fi

. pbs-env.sh

# Local goby locations
GOBY_DIR=${GOBY_DIR:-$SCRIPT_DIR/../..}

if [ ! -e ${GOBY_DIR}/goby.jar ]; then
    echo "goby.jar not found!"
    exit 2
fi

# The job name identifies the location for the scripts
JOB_TAG=$1
JOB_DIR=$(readlink -f .)/${JOB_TAG}
JOB_RESULTS_DIR=${JOB_RESULTS_DIR:-$(readlink -f .)/${JOB_TAG}-results}

if [ -e ${JOB_DIR} ]; then
    echo "The job output directory already exists"
    exit 3
fi

if [ ! -r ${REFERENCE} ]; then
    echo "WARNING: Reference ${REFERENCE} cannot be read"
fi

# Get the number of bytes in the reads file
if [ ! -r ${READS} ]; then
    echo "WARNING: Reads ${READS} file cannot be read"
    READS_SIZE=0
else
    READS_SIZE=`/usr/bin/stat --format=%s ${READS}`
fi

echo "Bundling job submission files"

# see if the job needs to split into a job array
if [ -z ${CHUNK_SIZE} ] || [ ${CHUNK_SIZE} -le 0 ] || [ ${CHUNK_SIZE} -ge ${READS_SIZE} ]; then
    echo "Alignment will not be split"
else
    NUMBER_OF_JOBS=$((${READS_SIZE} / ${CHUNK_SIZE} + 1))
    PBS_ARRAY_END_INDEX=$((${NUMBER_OF_JOBS} - 1))
    PBS_ARRAY_DIRECTIVE="#PBS -J0-${PBS_ARRAY_END_INDEX}"
    echo "Alignment will run as ${NUMBER_OF_JOBS} jobs"
fi

if [ ! -z ${PBS_STATUS_MAILTO} ]; then
    PBS_MAILTO_DIRECTIVE="#PBS -M ${PBS_STATUS_MAILTO}"
fi

# Use a basename that is based on name of the reads file
BASENAME=$(basename $READS .compact-reads)

# Create a somewhat unique tag but limit it to 15 characters (PBS name limit)
# ("stage" of job will be added - index, align, concat)
PBS_JOB_NAME=${JOB_TAG:0:10}

# Copy goby and submission scripts to the run directory
/bin/mkdir -p ${JOB_DIR}
/bin/cp ${GOBY_DIR}/goby.jar ${SCRIPT_DIR}/align.sh \
    ${SCRIPT_DIR}/index.sh ${SCRIPT_DIR}/align-concat.sh \
    ${SCRIPT_DIR}/index-align.sh ${SCRIPT_DIR}/concat.sh \
    ${SCRIPT_DIR}/index-align-concat.sh ${SCRIPT_DIR}/pbs-env.sh \
    ${JOB_DIR}

# Create job specific scripts from the template files
for FILE in goby-index.qsub goby-align.qsub goby-concat.qsub; do
    sed -e "s|%REFERENCE%|${REFERENCE}|" -e "s|%READS%|${READS}|" \
        -e "s|%PBS_QUEUE%|${PBS_QUEUE}|" -e "s|%PBS_MEMORY%|${PBS_MEMORY}|" \
        -e "s|%PBS_JVM_FLAGS%|${PBS_JVM_FLAGS}|" \
        -e "s|%PBS_ARRAY_DIRECTIVE%|${PBS_ARRAY_DIRECTIVE}|" \
        -e "s|%PBS_MAILTO_DIRECTIVE%|${PBS_MAILTO_DIRECTIVE}|" \
        -e "s|%REFERENCE_INDEX_NAME%|${REFERENCE_INDEX_NAME}|" \
        -e "s|%REFERENCE_INDEX_DIRECTORY%|${REFERENCE_INDEX_DIRECTORY}|" \
        -e "s|%ALIGNER%|${ALIGNER}|" -e "s|%COLORSPACE%|${COLORSPACE}|" \
        -e "s|%CHUNK_SIZE%|${CHUNK_SIZE}|" -e "s|%BASENAME%|${BASENAME}|" \
        -e "s|%PBS_JOB_NAME%|${PBS_JOB_NAME}|" \
        ${SCRIPT_DIR}/templates/${FILE} > ${JOB_DIR}/${FILE}
done

echo "Scripts were written to ${JOB_DIR}"
