#!/bin/bash
# Step 4: Convert BED to bedgraph and call peaks with SEACR.
# Reference: https://yezhengstat.github.io/CUTTag_tutorial/
#
#$ -cwd
#$ -o jobout.$JOB_ID.step4_peak_call
#$ -e joberr.$JOB_ID.step4_peak_call
#$ -j n
#$ -l h_rt=24:00:00,h_data=25G,highp
#$ -pe shared 1
#$ -M YOUR_EMAIL@institution.edu
#$ -m bea

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "$SCRIPT_DIR/config.sh"



run_loop() {
    f=$1; # feature
    b=$2; # biological replicate
    c=$3; # cell type

    histName="$f-$b-$c";

    subdir=${CUT_AND_TAG_DIR}/${f};
    bed_dir=$subdir/bed;
    bedgraph_dir=$subdir/bedgraph;
    seacr_dir=$subdir/peakCalling/SEACR;

    bed_fragments=$bed_dir/${histName}_bowtie2.fragments.bed;
    bedgraph=$bedgraph_dir/${histName}_bowtie2.fragments.normalized.bedgraph;
    
    mkdir -p $seacr_dir/log;
    mkdir -p $bedgraph_dir;

    cd $seacr_dir; # SEACR is gonna dump all the $password.auc and $password.auc.bed files here for some reason. not even a file name option, just "$password"???

    ### bed to bedgraph with genome coverages ###
    bedtools genomecov -bg -i $bed_fragments -g $GENOME_SIZES > $bedgraph;

    ### SEACR Peak Calling ###
    bash $SEACR $bedgraph 0.01 norm stringent $seacr_dir/${histName}_seacr_NoBg_T1.peaks &> $seacr_dir/log/${histName}_seacr_NoBg_T1.log;
    #mv $seacr_dir/${histName}_seacr_NoBg_T1.peaks.stringent.bed $seacr_dir/${histName}_seacr_NoBg_T1.peaks.stringent.no_blacklist.bed;
    bash $SEACR $bedgraph 0.1 norm stringent $seacr_dir/${histName}_seacr_NoBg_T10.peaks &> $seacr_dir/log/${histName}_seacr_NoBg_T10.log;
    #mv $seacr_dir/${histName}_seacr_NoBg_T10.peaks.stringent.bed $seacr_dir/${histName}_seacr_NoBg_T10.peaks.stringent.no_blacklist.bed;

    # remove blacklisted peaks
    bedtools intersect -a $seacr_dir/${histName}_seacr_NoBg_T1.peaks.stringent.bed -b $BED_BLACKLIST -v > $seacr_dir/${histName}_seacr_NoBg_T1.peaks.stringent.no_blacklist.bed;
    bedtools intersect -a $seacr_dir/${histName}_seacr_NoBg_T10.peaks.stringent.bed -b $BED_BLACKLIST -v > $seacr_dir/${histName}_seacr_NoBg_T10.peaks.stringent.no_blacklist.bed;
    
    # chr1-22, X, Y only
    awk '$1 ~ /^chr[0-9XY]/' $seacr_dir/${histName}_seacr_NoBg_T1.peaks.stringent.no_blacklist.bed > $seacr_dir/${histName}_seacr_NoBg_T1.peaks.stringent.no_blacklist.23chr_only.bed;
    awk '$1 ~ /^chr[0-9XY]/' $seacr_dir/${histName}_seacr_NoBg_T10.peaks.stringent.no_blacklist.bed > $seacr_dir/${histName}_seacr_NoBg_T10.peaks.stringent.no_blacklist.23chr_only.bed;

    ## sanity check
    #bash $SEACR $bedgraph 0.01 non stringent $seacr_dir/${histName}_seacr_NoBg_T1.non.peaks &> $seacr_dir/log/${histName}_seacr_NoBg_T1.non.log;
    #bash $SEACR $bedgraph 0.1 non stringent $seacr_dir/${histName}_seacr_NoBg_T10.non.peaks &> $seacr_dir/log/${histName}_seacr_NoBg_T10.non.log;

    ## sanity check 2
    #bash $SEACR $bedgraph 0.01 non relaxed $seacr_dir/${histName}_seacr_NoBg_T1.relaxed.peaks &> $seacr_dir/log/${histName}_seacr_NoBg_T1.relaxed.log;
    #bash $SEACR $bedgraph 0.1 non relaxed $seacr_dir/${histName}_seacr_NoBg_T10.relaxed.peaks &> $seacr_dir/log/${histName}_seacr_NoBg_T10.relaxed.log;
}


for f in $features; do
    for b in $bioReplicates; do
        for c in $cell_types; do
            run_loop $f $b $c;
        done
    done
done







echo "Step 4 (peak calling) completed."