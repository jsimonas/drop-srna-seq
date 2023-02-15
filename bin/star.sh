#!/bin/bash

set -euo pipefail

# run STAR
STAR \
    --genomeDir ${index} \
    --readFilesIn ${cdna_read} ${bc_read} \
    --soloCBwhitelist ${whitelist} \
    --runThreadN ${task.cpus} \
    --outFileNamePrefix ${prefix}_ \
    --readFilesCommand zcat \
    --runDirPerm All_RWX \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM GX GN \\
    --outSAMunmapped Within \
    --alignEndsType EndToEnd \
    --alignIntronMax 1 \
#    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMismatchNmax 1 \
    --outFilterMultimapScoreRange 0 \
    --outFilterMultimapNmax 50 \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outFilterMismatchNoverLmax 0.05
    --outFilterMatchNmin 16 \
#    --soloType CB_samTagOut \
    --soloType CB_UMI_Simple \
    --soloMultiMappers EM \
    --soloFeatures GeneFull \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 8 \
    --soloCBmatchWLtype 1MM