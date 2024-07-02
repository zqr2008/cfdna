### Required Arguments ###

bamFile=$1 # <*.mkdup.bam>

### Optional Arguments ###

outdir=$2
if [ ! -n "$outdir" ]; then
    outdir=`pwd`;
fi

LABEL=$3

MIN_MAPQ=$4
if [ ! -n "$MIN_MAPQ" ]; then
    MIN_MAPQ=30;
fi

MIN_MISMATCH=$5
if [ ! -n "$MIN_MISMATCH" ]; then
    MIN_MISMATCH=5;
fi

EXTEND_SIZE=$6
if [ ! -n "$EXTEND_SIZE" ]; then
    EXTEND_SIZE=10;
fi

### Preperation and Software ###

prefix=$(echo `basename $bamFile` | awk -F "." '{print $1}')
if [ ! -n "$LABEL" ]; then
    LABEL=$prefix
fi

bamtools=/path/bamtools
samtools=/path/samtools
perl=/usr/bin/perl
bindir=/path/
STEP1_1=${bindir}/bam_to_filtered_paired_sam_F12_V3.pl
STEP1_2=${bindir}/sam2bed_V2.pl
STEP2=${bindir}/bed_to_end_motif.pl
STEP3=${bindir}/bed_to_size_and_4mer_motif.short.long.peak.pl
STEP31=${bindir}/bed_to_size_and_preferend.v2.pl
STEP4=${bindir}/calcu_size_ratio_from_size_profile.pl
STEP5=${bindir}/calcu_motif_entropy_MDS.new.pl

check() {
    bamFile=$1
    checked=`$samtools quickcheck -v $bamFile && echo "Good" || echo "Fail"`
    if [ $checked = "Good" ]; then
        echo "PASS: $bamFile had been checked." > ${bamFile}.finish;
    elif [ $checked = "Fail" ]; then
        echo "FAIL: $bamFile failed the check.";
        exit
    else
        echo "ERROR: Checking function failed for $bamFile";
        exit
    fi
}

############
### Main ###
############

if [ -s "${outdir}/${prefix}.FILTER.finish" ]; then
    echo "[Exit] This Step had finished."
    exit
fi

if [ ! -s "$bamFile" ]; then
    echo "[Error] File Not Found: $bamFile does not exist."
    exit
fi

# Step 1-1. filter bam
#   1) remove duplicates
#   2) mapping quality >= MIN_MAPQ
#   3) only allow mismatch and mismatch (edit distance) <= MIN_MISMATCH
#   4) primary map, no multiple hits (by read name)
#
# Step 1-2 sam to fragment bed

if [ -e "${outdir}/${prefix}.fragment.bed.finish" ]; then
    echo "[Skip] Filter: ${outdir}/${prefix}.fragment.bed.finish exists."
else
    echo "[Filter] Start: `date`"

    $perl $STEP1_1 $bamFile $MIN_MAPQ $MIN_MISMATCH | \
    $perl $STEP1_2 - \
        > ${outdir}/${prefix}.fragment.bed && \
    touch ${outdir}/${prefix}.fragment.bed.finish

    echo "[Filter] End: `date`"
fi
# Step2. extract 10mer motif from fragment bed

if [ -e "${outdir}/${prefix}.motif.bed.finish" ]; then
    echo "[Skip] 10mer Motif: ${outdir}/${prefix}.motif.bed.finish exists."
elif [ ! -s "${outdir}/${prefix}.fragment.bed" ]; then
    echo "[Error] File Not Found: ${outdir}/${prefix}.fragment.bed does not exist."
    exit
else
    $perl $STEP2 $EXTEND_SIZE ${outdir}/${prefix}.fragment.bed \
        > ${outdir}/${prefix}.motif.bed && \
    rm -f ${outdir}/${prefix}.fragment.bed && \
    touch ${outdir}/${prefix}.motif.bed.finish
fi

# Step3. calcu size freq and motif freq from <*.motif.bed>

if [ -e "${outdir}/${prefix}.size.motif.finish" ]; then
    echo "[Skip] Size and 4mer Motif"
elif [ ! -s "${outdir}/${prefix}.motif.bed" ]; then
    echo "[Error] File Not Found: ${outdir}/${prefix}.motif.bed does not exist."
    exit
else
    $perl $STEP3 ${outdir}/${prefix}.motif.bed ${outdir}/${prefix} && \
    $perl $STEP31 ${outdir}/${prefix}.motif.bed ${outdir}/${prefix} && \	
    touch ${outdir}/${prefix}.size.motif.finish
fi

# Step4. combine size ratio together

if [ -e "${outdir}/${prefix}.size.finish" ]; then
    echo "[Skip] Size Ratio: ${outdir}/${prefix}.size.finish exists."
elif [ ! -s "${outdir}/${prefix}.size" ]; then
    echo "[Error] File Not Found: ${outdir}/${prefix}.size does not exist."
    exit
else
    $perl $STEP4 $prefix ${outdir}/${prefix}.size \
        > ${outdir}/${prefix}.size_ratio && \
    touch ${outdir}/${prefix}.size_ratio.finish
fi

# Step5. calculate motif diversity score (MDS).

if [ -e "${outdir}/${prefix}.4mer.motif.MDS.finish" ]; then
    echo "[Skip] Motif MDS: ${outdir}/${prefix}.4mer.motif.MDS.finish exists."
elif [ ! -s "${outdir}/${prefix}.4mer.motif.freq" ]; then
    echo "[Error] File Not Found: ${outdir}/${prefix}.4mer.motif.freq is required."
    exit
else
    $perl $STEP5 $LABEL ${prefix} ${outdir}/${prefix}.merge.4mer.motif.freq \
        > ${outdir}/${prefix}.merge.4mer.motif.MDS && \
    touch ${outdir}/${prefix}.4mer.motif.MDS.finish
fi

# Final check

if [ -s "${outdir}/${prefix}.filtered.bam.finish" ] && \
   [ -e "${outdir}/${prefix}.fragment.bed.finish" ] && \
   [ -e "${outdir}/${prefix}.motif.bed.finish" ] && \
   [ -e "${outdir}/${prefix}.size.motif.finish" ] && \
   [ -e "${outdir}/${prefix}.size_ratio.finish" ] && \
   [ -e "${outdir}/${prefix}.4mer.motif.MDS.finish" ]; then
    echo "Practice_makes_perfect" > ${outdir}/${prefix}.FILTER.finish
    rm ${outdir}/${prefix}.filtered.bam.finish ${outdir}/${prefix}.fragment.bed.finish ${outdir}/${prefix}.motif.bed.finish ${outdir}/${prefix}.size.motif.finish ${outdir}/${prefix}.size_ratio.finish ${outdir}/${prefix}.4mer.motif.MDS.finish
fi
