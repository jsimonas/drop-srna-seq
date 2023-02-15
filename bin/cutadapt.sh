#!/bin/bash

set -euo pipefail

# run cutadapt
cutadapt –u 4 –a A{8} -m 15 -j "$task.cpus" "${R2}.fastq" > "${R2}.fastq"