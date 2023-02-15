#!/bin/bash

set -euo pipefail

# run custom py script
convert_to_samplesheet.py --file "${sheet}" --out "standard_samplesheet.csv"
