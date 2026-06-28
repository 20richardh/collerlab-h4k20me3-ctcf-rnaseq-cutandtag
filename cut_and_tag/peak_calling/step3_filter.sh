#!/bin/bash
# Step 3: Filter alignments, convert to BED, and bin fragments.
# Reference: https://yezhengstat.github.io/CUTTag_tutorial/
#
#$ -cwd
#$ -o jobout.$JOB_ID.step3_filter
#$ -e joberr.$JOB_ID.step3_filter
#$ -j n
#$ -l h_rt=24:00:00,h_data=25G,highp
#$ -pe shared 1
#$ -M YOUR_EMAIL@institution.edu
#$ -m bea

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "$SCRIPT_DIR/config.sh"



### ALIGNMENT RESULTS FILTERING AND FILE FORMAT CONVERSION ###

run_loop() {
    f=$1; # feature
    b=$2; # biological replicate
    c=$3; # cell type
    histName="$f-$b-$c";
    
    bam_unsorted=$bam_dir/${histName}_bowtie2.bam;
    bam_sorted_rmDup=$bam_dir/${histName}_bowtie2.sorted.rmDup.bam;
    fragmentLen=$fragment_dir/${histName}_fragmentLen.txt;

    bam_sorted_rmDup_mappedOnly=$bam_dir/${histName}_bowtie2.sorted.rmDup.mapped_only.bam;
    bam_sorted_rmDup_mappedOnly_noBlacklist=$bam_dir/${histName}_bowtie2.sorted.rmDup.mapped_only.no_blacklist.bam;
    bam_sortedByName_rmDup_mappedOnly=$bam_dir/${histName}_bowtie2.sortedByName.rmDup.mapped_only.bam;
    bed_raw=$bed_dir/${histName}_bowtie2.sorted.rmDup.mapped_only.raw.bed;
    bed_clean=$bed_dir/${histName}_bowtie2.sorted.rmDup.mapped_only.clean.bed;
    bed_clean_chr=$bed_dir/${histName}_bowtie2.sorted.rmDup.mapped_only.clean.23chr_only.bed;
    bed_fragments=$bed_dir/${histName}_bowtie2.fragments.bed;
    bed_fragments_chr=$bed_dir/${histName}_bowtie2.fragments.23chr_only.bed;
    bin_len_bed=$bed_dir/${histName}_bowtie2.fragmentsCount.bin$BIN_LEN.bed;
    bin_len_bed_chr=$bed_dir/${histName}_bowtie2.fragmentsCount.bin$BIN_LEN.23chr_only.bed;

    ## Extract the 9th column from the alignment sam file which is the fragment length
    samtools view -F 0x04 $bam_sorted_rmDup | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > $fragmentLen;
    ## Filter and keep the mapped read pairs
    samtools view -bS -F 0x04 $bam_sorted_rmDup > $bam_sorted_rmDup_mappedOnly;

    ### Sort by name, due to bamtobed requirements
    # NOTE 1: made a file for bam reads without blacklist regions. Could implement for peak calling but doesn't change much---------------------------------------------
    # NOTE 2: samtools sort -n is used to sort by read name, which is required for bedtools bamtobed to work properly.---------------------------------------------
    bedtools intersect -abam $bam_sorted_rmDup -b $BED_BLACKLIST -v | samtools view -bS -F 0x04 - > $bam_sorted_rmDup_mappedOnly_noBlacklist;
    samtools index $bam_sorted_rmDup_mappedOnly;
    samtools index $bam_sorted_rmDup_mappedOnly_noBlacklist;
    samtools sort -n $bam_sorted_rmDup_mappedOnly -o $bam_sortedByName_rmDup_mappedOnly;
    # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    ## Convert into bed file format
    bedtools bamtobed -i $bam_sortedByName_rmDup_mappedOnly -bedpe > $bed_raw;
    ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    ## NOTE 3: ALSO use only chr1-22 and chrX and chrY (BJ fibroblasts, all male) -- Richard--------------------------------------------------------------------------------------------
    awk '$1==$4 && $6-$2 < 1000' $bed_raw > $bed_clean;
    awk '$1 ~ /^chr[0-9XY]/ {print $0}' $bed_clean > $bed_clean_chr;
    # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Only extract the fragment related columns
    cut -f 1,2,6 $bed_clean | sort -k1,1 -k2,2n -k3,3n > $bed_fragments;
    cut -f 1,2,6 $bed_clean_chr | sort -k1,1 -k2,2n -k3,3n > $bed_fragments_chr;

    ## We use the mid point of each fragment to infer which [BIN_LEN] bp bins does this fragment belong to.
    awk -v w=$BIN_LEN '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $bed_fragments | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  > $bin_len_bed;
    awk -v w=$BIN_LEN '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $bed_fragments_chr | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  > $bin_len_bed_chr;

    ## do 500bp as well
    awk -v w=500 '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $bed_fragments | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n > $bed_dir/${histName}_bowtie2.fragmentsCount.bin500.bed;
    awk -v w=500 '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $bed_fragments_chr | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n > $bed_dir/${histName}_bowtie2.fragmentsCount.bin500.23chr_only.bed;

    ## Clean up
    #rm $bed_raw;
}

for f in $features; do
    subdir=${CUT_AND_TAG_DIR}/${f};
    bam_dir=$subdir/bam;
    fragment_dir=$bam_dir/fragmentLen;
    bed_dir=$subdir/bed;
    
    mkdir -p $fragment_dir;
    mkdir -p $bed_dir;
    for b in $bioReplicates; do
        for c in $cell_types; do
            run_loop $f $b $c;
        done
    done
done



echo "Step 3 (filter and convert) completed."