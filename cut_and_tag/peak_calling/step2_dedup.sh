#!/bin/bash
# Step 2: Sort BAM files and remove duplicates with Picard.
# Reference: https://yezhengstat.github.io/CUTTag_tutorial/
#
#$ -cwd
#$ -o jobout.$JOB_ID.step2_dedup
#$ -e joberr.$JOB_ID.step2_dedup
#$ -j n
#$ -l h_rt=24:00:00,h_data=25G,highp
#$ -pe shared 4
#$ -M YOUR_EMAIL@institution.edu
#$ -m bea

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "$SCRIPT_DIR/config.sh"

picardCMD="java -jar $PICARD"



### Removing Duplicates using PICARD ###
for f in $features; do
    subdir=${CUT_AND_TAG_DIR}/${f};
    fastq_dir=$subdir/fastq;
    bam_dir=$subdir/bam;
    bed_dir=$subdir/bed;
    bedgraph_dir=$subdir/bedgraph;

    mkdir -p $bam_dir/picard_summary;

    for b in $bioReplicates; do
        for c in $cell_types; do
            histName="$f-$b-$c";
            bam_unsorted=${bam_dir}/${histName}_bowtie2.bam;
            bam_sorted=${bam_dir}/${histName}_bowtie2.sorted.bam;
            bam_sorted_dupMark=${bam_dir}/${histName}_bowtie2.sorted.dupMarked.bam;
            metrics_dupMark=$bam_dir/picard_summary/${histName}_picard.dupMark.txt;
            bam_sorted_rmDup=${bam_dir}/${histName}_bowtie2.sorted.rmDup.bam;
            metrics_rmDup=$bam_dir/picard_summary/${histName}_picard.rmDup.txt;

            # PICARD sort by coordinates
            $picardCMD SortSam I=$bam_unsorted O=$bam_sorted SORT_ORDER=coordinate;
            # PICARD mark and remove duplicates
            $picardCMD MarkDuplicates I=$bam_sorted O=$bam_sorted_dupMark METRICS_FILE=$metrics_dupMark;
            $picardCMD MarkDuplicates I=$bam_sorted O=$bam_sorted_rmDup  REMOVE_DUPLICATES=true METRICS_FILE=$metrics_rmDup;
        done
    done
done

echo "Step 2 (deduplication) completed."