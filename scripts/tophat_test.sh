#!/usr/bin/env bash
# Tophat SL mapping test
CONDITION="mismatches-0_minlength-10"
BASE_OUTPUTDIR=build/tophat/${CONDITION}
GENOME_DIR=/cbcb/lab/nelsayed/ref_data/tcruzi_clbrener/genome/tc_esmer
REF=${GENOME_DIR}/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome
TOPHAT_ARGS="--num-threads 8 --max-multihits 1 --mate-inner-dist 170"

INPUT_DIR=build/01-filtered_reads/${CONDITION}
R1_INPUT=${INPUT_DIR}/HPGL0258_R1_combined.filtered.fastq.gz
R2_INPUT=${INPUT_DIR}/HPGL0258_R2_combined.filtered.fastq.gz
R2_INPUT_SE=${INPUT_DIR}/HPGL0258_R2_combined.filtered.fastq.gz
R1_INPUT_UNTRIMMED=${INPUT_DIR}/HPGL0258_R1_combined.filtered_complete.fastq.gz
R2_INPUT_UNTRIMMED=${INPUT_DIR}/HPGL0258_R2_combined.filtered_SE_complete.fastq.gz

#######################################
# SL False hits (R1)
#######################################
echo "Finding false SL hits in matched R1 reads"
OUTPUTDIR=${BASE_OUTPUTDIR}/R1_false_hits
mkdir -p ${OUTPUTDIR}

TOPHAT_CMD="/usr/bin/tophat $TOPHAT_ARGS  \
                   --no-mixed             \
                   -o ${OUTPUTDIR}        \
                   $REF                   \
                   $R1_INPUT_UNTRIMMED    \
                   $R2_INPUT"
echo $TOPHAT_CMD
eval $TOPHAT_CMD

# Sort and index output
samtools sort -@ 12 ${OUTPUTDIR}/accepted_hits.bam \
                    ${OUTPUTDIR}/accepted_hits_sorted
samtools index ${OUTPUTDIR}/accepted_hits_sorted.bam

#######################################
# SL False hits (R2)
#######################################
echo "Finding false SL hits in matched R2 reads"
OUTPUTDIR=${BASE_OUTPUTDIR}/R2_false_hits
mkdir -p ${OUTPUTDIR}

TOPHAT_CMD="/usr/bin/tophat $TOPHAT_ARGS  \
                   -o ${OUTPUTDIR} \
                   $REF          \
                   $R2_INPUT_UNTRIMMED"
echo $TOPHAT_CMD
eval $TOPHAT_CMD

# Sort and index output
samtools sort -@ 12 ${OUTPUTDIR}/accepted_hits.bam \
                    ${OUTPUTDIR}/accepted_hits_sorted
samtools index ${OUTPUTDIR}/accepted_hits_sorted.bam

#######################################
# SL found in R1 (R1 + R2)
#######################################
echo "Processing R1 SL matches"
OUTPUTDIR=${BASE_OUTPUTDIR}/R1_matches
mkdir -p ${OUTPUTDIR}

TOPHAT_CMD="/usr/bin/tophat $TOPHAT_ARGS  \
                   -o ${OUTPUTDIR} \
                   $REF          \
                   $R1_INPUT     \
                   $R2_INPUT"
echo $TOPHAT_CMD
eval $TOPHAT_CMD

# Sort and index output
samtools sort -@ 12 ${OUTPUTDIR}/accepted_hits.bam \
                    ${OUTPUTDIR}/accepted_hits_sorted
samtools index ${OUTPUTDIR}/accepted_hits_sorted.bam

#######################################
# SL found in R2 (R2 only)
#######################################
echo "Processing R2 SL matches"
OUTPUTDIR=${BASE_OUTPUTDIR}/R2_matches
mkdir -p ${OUTPUTDIR}

TOPHAT_CMD="/usr/bin/tophat $TOPHAT_ARGS  \
                   -o ${OUTPUTDIR} \
                   $REF          \
                   $R2_INPUT_SE"
echo $TOPHAT_CMD
eval $TOPHAT_CMD

# Sort and index output
samtools sort -@ 12 ${OUTPUTDIR}/accepted_hits.bam \
                    ${OUTPUTDIR}/accepted_hits_sorted
samtools index ${OUTPUTDIR}/accepted_hits_sorted.bam

