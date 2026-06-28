#!/bin/bash
# =============================================================================
# USER CONFIGURATION — set PROJ_PATH and EMAIL before running the pipeline.
#
# Expected directory layout under PROJ_PATH:
#   data/genome/              reference genome + bowtie2 index
#   data/cut_and_tag/<mark>/  fastq/, bam/, bed/, bedgraph/ per histone mark
#
# Scheduler: scripts use SGE (#$ directives) as used on UCLA Hoffman2.
# Adapt #$ directives and module load commands for your HPC environment.
# =============================================================================

PROJ_PATH=/path/to/your/project
EMAIL=YOUR_EMAIL@institution.edu

# --- Derived paths (adjust only if your layout differs) ----------------------
DATA_DIR=$PROJ_PATH/data

# Module loading — comment out if not needed, or adjust module names/versions for HPC environment
. /u/local/Modules/default/init/modules.sh
module load bowtie2/2.4.2
module load fastqc/0.11.9
module load picard_tools
module load samtools/1.15
module load bedtools/2.30.0
module load R

# Conda environment providing deeptools
source ~/anaconda3/cut_and_tag/bin/activate

# --- Experiment metadata -----------------------------------------------------
features="H4K20me3 CTCF"
bioReplicates="12-3 12-4"
cell_types="P CI SS"

# --- Auto-detected repo paths ------------------------------------------------
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
REPO_DIR=$(realpath "$SCRIPT_DIR/../..")
BED_BLACKLIST=$REPO_DIR/data/reference/GRCh38_unified_blacklist.bed
SEACR=$SCRIPT_DIR/SEACR_1.3.sh

# --- Genome files (large external data under PROJ_PATH) ----------------------
GENOME_DIR=$DATA_DIR/genome
GENOME_NAME=GRCh38.primary_assembly.genome
GENOME=$GENOME_DIR/$GENOME_NAME.fa
GENOME_GTF=$GENOME_DIR/gencode.v29.primary_assembly.annotation.gtf
GENOME_SIZES=$GENOME_DIR/$GENOME_NAME.sizes
BOWTIE2_INDEX=$GENOME_DIR/bowtie2/hg38_index/hg38_index

# --- CUT&Tag directories -----------------------------------------------------
CUT_AND_TAG_DIR=$DATA_DIR/cut_and_tag

# --- Pipeline parameters -----------------------------------------------------
MIN_QUALITY_SCORE=2
BIN_LEN=200
