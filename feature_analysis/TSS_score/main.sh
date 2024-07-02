#!/bin/bash
# Author: Li lingguo
# Usage: sh main.sh <bam_file> <tss_bed> <ref> <depth> <outdir> <chr> 


# Required Arguments
bam_file=$1
tss_bed=$2
ref=$3
depth=$4
CHR=$6

# Optional Arguments
outdir=$5

if [ ! -n "${outdir}" ]; then
    outdir=$(pwd)
fi
if [ ! -d ${outdir} ]; then
    mkdir -p ${outdir}
fi


# Preparation
pre_prefix=$(echo $(basename $bam_file) | awk -F "." '{print $1}')
prefix=$(echo "${pre_prefix}.chr${CHR}")

# Environment and Softwares
samtools=/path/samtools
python3=/path/python3

get_score=/TSS_path/mpileup_to_tss_score.py
get_coverage=/TSS_path/mpileup_to_coverage.py
get_binDep=/TSS_path/coverage_to_bin_depth.py
get_TSS=/TSS_path/binDep_to_Tss.py

### Main ###

##step1 mpileup
if [ -s ${outdir}/${prefix}.finish ]; then
    echo "[Skip] mpileup: ${outdir}/${prefix}.finish exists."
else
    ${samtools} mpileup -B -f ${ref} -l ${tss_bed} ${bam_file} > ${outdir}/${prefix}.mpileup
    echo "Finish" > ${outdir}/${prefix}.finish
fi

## step2 mpileup to core TSS_score
if [ ! -s "${outdir}/${prefix}.finish" ]; then
    echo "[Error] Score: ${outdir}/${prefix}.finish does not exist."
    exit
elif [ -s "${outdir}/${prefix}.score" ]; then
    echo "[Skip] Score：${outdir}/${prefix}.score exists."
else
    ${python3} ${get_score} ${tss_bed} ${outdir}/${prefix}.mpileup ${depth} ${outdir}/${prefix}.score 1000
fi

## step3 mpileup to tss extend coverage
if [ ! -s ${outdir}/${prefix}.finish ]; then
    echo "[Error] coverage: ${outdir}/${prefix}.finish does not exist."
    exit
elif [ -s "${outdir}/${prefix}.coverage" ]; then
    echo "[Skip] Coverage: ${outdir}/${prefix}.coverage exits."
else
    ${python3} ${get_coverage} ${tss_bed} ${outdir}/${prefix}.mpileup ${depth} ${outdir}/${prefix}.coverage
fi

## step4 tss coverage to 10_bin_Depth
if [ ! -s ${outdir}/${prefix}.finish ]; then
    echo "[Error] binDep: ${outdir}/${prefix}.finish does not exist."
    exit
elif [ ! -s ${outdir}/${prefix}.coverage ]; then
    echo "[Error] binDep: ${outdir}/${prefix}.coverage does not exist."
    exit
elif [ -s "${outdir}/${prefix}.binDep" ]; then
    echo "[Skip] binDep: ${outdir}/${prefix}.binDep exits."
else
    ${python3} ${get_binDep}  ${outdir}/${prefix}.coverage  ${outdir}/${prefix}.binDep
fi

## step5 10_bin_Depth to 1_TSS_score
if [ ! -s ${outdir}/${prefix}.finish ]; then
    echo "[Error] TSS: ${outdir}/${prefix}.finish does not exist."
    exit
elif [ ! -s ${outdir}/${prefix}.binDep ]; then
    echo "[Error] TSS: ${outdir}/${prefix}.binDep does not exist."
    exit
elif [ -s "${outdir}/${prefix}.TSS_score" ]; then
    echo "[Skip] TSS: ${outdir}/${prefix}.TSS_score exits."
else
    ${python3} ${get_TSS}  ${outdir}/${prefix}.binDep  ${outdir}/${prefix}.TSS_score
fi
