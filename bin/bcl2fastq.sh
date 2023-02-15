#!/bin/bash

set -euo pipefail

# run bcl2fastq
 bcl2fastq \
    --runfolder-dir ${runDir} \
    --output-dir . \
    --sample-sheet ${sheet} \
    --mask-short-adapter-reads 0 \
    --minimum-trimmed-read-length 0 \
    --use-bases-mask 'y*,I*,y*,y*' \
    --no-lane-splitting \
    --create-fastq-for-index-reads \
    --processing-threads $task.cpus