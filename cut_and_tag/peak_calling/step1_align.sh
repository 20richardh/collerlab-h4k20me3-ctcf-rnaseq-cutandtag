#!/bin/bash
# Step 1: Align FASTQ files to the genome with Bowtie2.
# Reference: https://yezhengstat.github.io/CUTTag_tutorial/
#
# SGE scheduler directives (UCLA Hoffman2). Adapt for your HPC environment.
#$ -cwd
#$ -o jobout.$JOB_ID.step1_align
#$ -e joberr.$JOB_ID.step1_align
#$ -j n
#$ -l h_rt=24:00:00,h_data=25G,highp
#$ -pe shared 4
#$ -M YOUR_EMAIL@institution.edu
#$ -m bea

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "$SCRIPT_DIR/config.sh"

# Before running, create hg38_index folder (see config.sh BOWTIE2_INDEX).
ref=$BOWTIE2_INDEX

cores=4


merge_fastq() {
    filenames=$1;
    outfile=$2;
    cat $filenames > $outfile;
}

for f in $features; do
    subdir=${CUT_AND_TAG_DIR}/${f}
    fastq_dir=$subdir/fastq
    bam_dir=$subdir/bam
    bed_dir=$subdir/bed
    bedgraph_dir=$subdir/bedgraph

    mkdir -p $bam_dir/bowtie2_summary
    mkdir -p $bed_dir
    mkdir -p $bedgraph_dir

    for b in $bioReplicates; do
        for c in $cell_types; do
            histName="$f-$b-$c";

            files_1=$(ls ${fastq_dir}/${b}_${c}_${f}_S*_R1_*.fastq.gz);
            files_2=$(ls ${fastq_dir}/${b}_${c}_${f}_S*_R2_*.fastq.gz);

            fastq_merged_1=${fastq_dir}/${b}_${c}_${f}_R1.merged.fastq.gz;
            fastq_merged_2=${fastq_dir}/${b}_${c}_${f}_R2.merged.fastq.gz;

            merge_fastq $files_1 $fastq_merged_1;
            merge_fastq $files_2 $fastq_merged_2;

            temp_sam=${bam_dir}/${histName}_bowtie2.sam;
            bam=${bam_dir}/${histName}_bowtie2.bam;

            bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 $fastq_merged_1 -2 $fastq_merged_2 -S $temp_sam &> ${bam_dir}/bowtie2_summary/${histName}_bowtie2.txt;
            samtools view -@ $cores -bS -o $bam $temp_sam;
            rm $temp_sam;
        done
    done
done

echo "Step 1 (alignment) completed."