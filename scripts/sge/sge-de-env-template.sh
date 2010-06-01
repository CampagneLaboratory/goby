#
# Environment variables used for SGE alignment jobs
#

#
# Email address to notify upon job completion
#
SGE_STATUS_MAILTO=

 #
 # Email address to notify upon job completion
 #
SGE_STATUS_MAILTO=fac2003@gmail.com

 #
 # Memory size for SGE

 #
SGE_MEMORY=5g

 #
 # JVM flags (memory should match SGE_MEMORY)
 #
SGE_JVM_FLAGS=-Xmx16g


alignment_basenames="\
      DLTTEJH-Bullard-HBR-SRR037439.entries  \
      HHNVSNR-Bullard-HBR-SRR037440.entries  \
      RRTFBOP-Bullard-HBR-SRR037443.entries  \
      ZOBXCNB-Bullard-HBR-SRR037441.entries  \
      DOWTGPI-Bullard-HBR-SRR037444.entries  \
      WKQCRQC-Bullard-HBR-SRR037442.entries  \
      LHFEDQE-solid-HBR.entries"

group1=HBR-Bullard-ILM
group2=HBR-SOLID
group1_basenames=DLTTEJH-Bullard-HBR-SRR037439.entries,DOWTGPI-Bullard-HBR-SRR037444.entries,HHNVSNR-Bullard-HBR-SRR037440.entries,RRTFBOP-Bullard-HBR-SRR037443.entries,WKQCRQC-Bullard-HBR-SRR037442.entries,ZOBXCNB-Bullard-HBR-SRR037441.entries
group2_basenames=LHFEDQE-solid-HBR
annotation_file=exon-annotations-NCBI36.54.tsv
use_weights=gc
adjust_gc_bias_boolean=true