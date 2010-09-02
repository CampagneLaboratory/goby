#!/bin/bash -x
#
# This script will help you
# 1. Index a reference for a Transcript alignment
# 2. Align a reads file against that reference
#

##
##  Edit the following items...
##
TAG=MeaningfulTag
COMPACT_READS_FILE_TO_ALIGN=your_reads_file.compact-reads
REFERENCE_FASTA_GZ=Homo_sapiens.GRCh37.57.cdna.all.fa.gz
VERSION=GRCh37.57
ORGANISM=homo_sapiens
JVM_FLAGS=-Xmx8g

#
# The aligner to use
#
ALIGNER=bwa
ALIGNER_OPTIONS=
QUALITY_FILTER_PARAMETERS="threshold=0.05"
AMBIGUITY_THRESHOLD=2

#
# When splitting the input READS file, how big to make each chunk. In bytes.
#
READS_SPLIT_CHUNK_SIZE=100000000

#
# Should be basespace or colorspace... depends on the reads you will be aligning
#
SPACE=basespace


############################################################
############################################################
##                                                        ##
##  You shouldn't have to change anything below here.     ##
##                                                        ##
############################################################
############################################################

function removeDir {
    if [ -d $1 ]; then
        rm -rf $1
    fi
}

function calculate_PAD_FORMAT {
    _NUMBER=$1
    _NUMBER_TO_PAD=0
    while [ ${_NUMBER} -gt 10 ]; do
        _NUMBER_TO_PAD=$(( _NUMBER_TO_PAD + 1 ))
        _NUMBER=$(( _NUMBER / 10 ))
    done
    if [ ${_NUMBER_TO_PAD} -eq 0 ]; then
        PAD_FORMAT=%d
    else
        PAD_FORMAT=%0$(( _NUMBER_TO_PAD + 1 ))d
    fi
}

function cleanup {
    # Cleanup failed previous attempts for this tag
    removeDir ${SPLIT_INDEX_DIR}
    removeDir ${SPLIT_ALIGN_DIR}
}

function init {
    echo "########################################################"
    echo "###  init"
    echo "########################################################"
    SPACE_PARAM=""
    if [ $SPACE = "colorspace" ]; then
        SPACE_PARAM="--color-space"
    fi
    
    INDEX_PREFIX=index
    CURRENT_DIR=`pwd`
    OUTPUT_BASE_DIR=${CURRENT_DIR}/reference-db
    SPLIT_WORK_DIR=${CURRENT_DIR}/work-${TAG}
    SPLIT_INDEX_DIR=${SPLIT_WORK_DIR}-index
    SPLIT_ALIGN_DIR=${SPLIT_WORK_DIR}-align
    REFERENCE_FASTA_GZ=${CURRENT_DIR}/${REFERENCE_FASTA_GZ}
    GOBY_JAR=${CURRENT_DIR}/goby.jar
    LOG4J_CONFIG=${CURRENT_DIR}/log4j.properties
    GOBY_CONFIG=${CURRENT_DIR}/goby.properties
    
    #
    # Rename version since this is a transcript index
    #
    VERSION=Transcript-${VERSION}
    
    #
    # Directory where we store the compact reference and toplevel-ids file.
    #
    INPUT_COMPACT_REF_DIR=${OUTPUT_BASE_DIR}/${VERSION}/${ORGANISM}/reference
    INDEXED_REF_DIR=${OUTPUT_BASE_DIR}/${VERSION}/${ORGANISM}/${SPACE}/${ALIGNER}

    # How many parts to split the READS file into
    READS_FILE_SIZE=$(stat -c%s "${COMPACT_READS_FILE_TO_ALIGN}")
    NUMBER_OF_ALIGN_PARTS=$((READS_FILE_SIZE / READS_SPLIT_CHUNK_SIZE))
    if [ ${READS_FILE_SIZE} -gt $(( NUMBER_OF_ALIGN_PARTS * READS_SPLIT_CHUNK_SIZE ))  ]; then
         NUMBER_OF_ALIGN_PARTS=$(( NUMBER_OF_ALIGN_PARTS + 1 ))
    fi

    PRE_RESULTS_DIR=${SPLIT_ALIGN_DIR}/pre-results
}

function indexReference {
    echo "########################################################"
    echo "###  indexReference"
    echo "########################################################"
    #
    # See if it is likely that the index has already been built. If so,
    # don't rebuild it.
    #
    if [ ! -f "${INPUT_COMPACT_REF_DIR}/gene-transcript-map.txt" ]; then
    
        removeDir ${INPUT_COMPACT_REF_DIR}
        removeDir ${INDEXED_REF_DIR}
    
        mkdir -p ${INPUT_COMPACT_REF_DIR}/
    
        #
        # Split the fasta by transcript
        #
        mkdir -p ${SPLIT_INDEX_DIR}
    
        java ${JVM_FLAGS} \
            -Dlog4j.configuration=file:${LOG4J_CONFIG} \
            -Dgoby.configuration=file:${GOBY_CONFIG} \
            -jar ${GOBY_JAR} --mode split-transcripts --input ${REFERENCE_FASTA_GZ} \
            --output ${SPLIT_INDEX_DIR}/${INDEX_PREFIX}
        RETURN_STATUS=$?
        if [ ! $RETURN_STATUS -eq 0 ];then
            echo "Split reference by transcript failed"
            exit
        fi
        mv ${SPLIT_INDEX_DIR}/${INDEX_PREFIX}.config ${INPUT_COMPACT_REF_DIR}/gene-transcript-map.txt
    
        cd ${SPLIT_INDEX_DIR}
        java ${JVM_FLAGS} \
            -Dlog4j.configuration=file:${LOG4J_CONFIG} \
            -Dgoby.configuration=file:${GOBY_CONFIG} \
            -jar ${GOBY_JAR} --mode fasta-to-compact -n 1 --include-identifiers *.fa.gz
        RETURN_STATUS=$?
        if [ ! $RETURN_STATUS -eq 0 ];then
            echo "Conversion of plit reference by transcript files to .compact-reads files failed."
            exit
        fi
        mv *.compact-reads ${INPUT_COMPACT_REF_DIR}
        rm *.fa.gz
        cd ${CURRENT_DIR}
    
        cp ${CURRENT_DIR}/goby.properties ${SPLIT_INDEX_DIR}
    
        for SUB_COMPACT_READS in ${INPUT_COMPACT_REF_DIR}/*.compact-reads
        do
            SUB_INDEX_PREFIX=`basename ${SUB_COMPACT_READS}`
            SUB_INDEX_PREFIX=${SUB_INDEX_PREFIX%\.*}T
            echo "Creating index ${ALIGNER} / ${SPACE} / ${SUB_INDEX_PREFIX} ..."
            cd ${SPLIT_INDEX_DIR}
            java ${JVM_FLAGS} \
                -Dlog4j.configuration=file:${LOG4J_CONFIG} \
                -Dgoby.configuration=file:${GOBY_CONFIG} \
                -jar ${GOBY_JAR} --mode align --index --reference ${SUB_COMPACT_READS} ${SPACE_PARAM} \
                --aligner ${ALIGNER} --database-name \
                ${INDEXED_REF_DIR}/${SUB_INDEX_PREFIX}
            RETURN_STATUS=$?
            if [ ! $RETURN_STATUS -eq 0 ];then
                echo "Indexing of transcript reference failed for ${SUB_COMPACT_READS}"
                exit
            fi
            cd ${CURRENT_DIR}
        done
    fi
}

#
# Split the reads
# Figure out how many transcript parts there were
#
function splitReads {
    echo "########################################################"
    echo "###  splitReads"
    echo "########################################################"
    NUMBER_OF_TRANSCRIPT_PARTS=`ls ${INPUT_COMPACT_REF_DIR}/${INDEX_PREFIX}.*.compact-reads | wc -l`

    #
    # Split the reads
    #  
    mkdir -p ${SPLIT_ALIGN_DIR}
    if [ ${NUMBER_OF_ALIGN_PARTS} -eq 1 ]; then
        cp ${COMPACT_READS_FILE_TO_ALIGN} ${SPLIT_ALIGN_DIR}/0.compact-reads
    else
        for ((i = 0; i < ${NUMBER_OF_ALIGN_PARTS}; i++));
        do
            READS_FILE=${SPLIT_ALIGN_DIR}/${i}.compact-reads
    
            START_POSITION=$(( i *  READS_SPLIT_CHUNK_SIZE ))
            END_POSITION=$(( START_POSITION +  READS_SPLIT_CHUNK_SIZE - 1 ))
    
            java ${JVM_FLAGS} \
                -Dlog4j.configuration=file:${LOG4J_CONFIG} \
                -Dgoby.configuration=file:${GOBY_CONFIG} \
                -jar ${GOBY_JAR} --mode reformat-compact-reads --output ${READS_FILE} \
                --start-position ${START_POSITION} --end-position ${END_POSITION} ${COMPACT_READS_FILE_TO_ALIGN}
            RETURN_STATUS=$?
            if [ ! $RETURN_STATUS -eq 0 ];then
                echo "Spliting the input reads file failed."
                exit
            fi
        done
    fi
}

#
# Align the split reads
#
function alignSplitReads {
    echo "########################################################"
    echo "###  alignSplitReads"
    echo "########################################################"
    calculate_PAD_FORMAT ${NUMBER_OF_TRANSCRIPT_PARTS}
    for ((READS_NUM = 0; READS_NUM < NUMBER_OF_ALIGN_PARTS; READS_NUM++));
    do
        READS_FILE=${SPLIT_ALIGN_DIR}/${READS_NUM}.compact-reads
        for ((TRANSCRIPT_NUMBER = 0; TRANSCRIPT_NUMBER < NUMBER_OF_TRANSCRIPT_PARTS; TRANSCRIPT_NUMBER++));
        do
            REF_INDEX_PREFIX=`printf ${INDEX_PREFIX}.${PAD_FORMAT} ${TRANSCRIPT_NUMBER}`
            REFERENCE_INDEX_FILE=${INPUT_COMPACT_REF_DIR}/${REF_INDEX_PREFIX}.compact-reads        
            java ${JVM_FLAGS} \
                -Dlog4j.configuration=file:${LOG4J_CONFIG} \
                -Dgoby.configuration=file:${GOBY_CONFIG} \
                -jar ${GOBY_JAR} \
                --mode align \
                --reference ${REFERENCE_INDEX_FILE} \
                --aligner ${ALIGNER} ${SPACE_PARAM} --search \
                --ambiguity-threshold ${AMBIGUITY_THRESHOLD} --quality-filter-parameters "${QUALITY_FILTER_PARAMETERS}" \
                --database-name ${REF_INDEX_PREFIX}T --database-directory ${INDEXED_REF_DIR} \
                ${ALIGNER_OPTIONS} --reads ${READS_FILE} --basename ${TAG}
             RETURN_STATUS=$?
             rm ${READS_NUM}.fasta ${READS_NUM}.fastq
             if [ ! $RETURN_STATUS -eq 0 ];then
                 echo "Alignment failed"
                 exit
             fi
    
             # Alignment completed for this part
             RESULT_DIR=${PRE_RESULTS_DIR}/transcripts-${READS_NUM}-${TRANSCRIPT_NUMBER}
             /bin/mkdir -p ${RESULT_DIR}
             /bin/mv *.entries *.header *.stats *.tmh ${RESULT_DIR}
        done
    done
}

#
# Merge the results
#
function merge {
    echo "########################################################"
    echo "###  merge"
    echo "########################################################"
    GENE_TRANSCRIPT_MAP_FILE=${INPUT_COMPACT_REF_DIR}/gene-transcript-map.txt
    for ((READS_NUM = 0; READS_NUM < NUMBER_OF_ALIGN_PARTS; READS_NUM++));
    do
        if [ ${NUMBER_OF_ALIGN_PARTS} -eq 1 ]; then
            OUTPUT_DIR=${CURRENT_DIR}/results-${TAG}
        else
            OUTPUT_DIR=${PRE_RESULTS_DIR}/${TAG}-${READS_NUM}
        fi
        mkdir -p ${OUTPUT_DIR}
        java ${JVM_FLAGS} \
            -Dlog4j.configuration=file:${LOG4J_CONFIG} \
            -Dgoby.configuration=file:${GOBY_CONFIG} \
            -jar ${GOBY_JAR} \
            --mode merge-compact-alignments \
            --gene-transcript-map-file ${GENE_TRANSCRIPT_MAP_FILE} \
            --output ${OUTPUT_PREFIX}/${TAG} \
            ${PRE_RESULTS_DIR}/transcripts-${READS_NUM}-*/*.entries
         RETURN_STATUS=$?
         if [ ! $RETURN_STATUS -eq 0 ];then
             echo "Merge failed"
             exit
         fi
    done
}

#
# Concat
#
function concat {
    echo "########################################################"
    echo "###  concat"
    echo "########################################################"
    if [ ${NUMBER_OF_ALIGN_PARTS} -gt 1 ];then
        RESULT_DIR=${CURRENT_DIR}/results-${TAG}
        mkdir -p ${RESULT_DIR}
        java ${JVM_FLAGS} \
            -Dlog4j.configuration=file:${LOG4J_CONFIG} \
            -Dgoby.configuration=file:${GOBY_CONFIG} \
            -jar ${GOBY_JAR} \
            --mode concatenate-alignments --adjust-query-indices false \
            --output ${RESULT_DIR}/${TAG} ${PRE_RESULTS_DIR}/${TAG}-*/*.entries
        RETURN_STATUS=$?
        if [ ! $RETURN_STATUS -eq 0 ];then
            echo "Concat failed"
            exit
        fi
    
    fi
}

#
# Execute the transcript alignment
#
cleanup
init
indexReference
splitReads
alignSplitReads
merge
concat
cleanup
