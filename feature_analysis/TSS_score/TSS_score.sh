## old useage
## sh step00.sh <input.list> <outdir> <main.sh> <genome.ref> <wgs.bed> <depth> <0 or 01-24>

## useage
## sh step00.sh <input.list> <outdir> 

list=$1	# <group:drname> <sampleID:sub dirname> <bam>
path=$(echo -e "$(dirname $2)/$(basename $2)")


#-----------------------------------------------------------------
shell=/TSS_path/main.sh
ref=Homo_sapiens_assembly38.fa
bed=/TSS_path/hg38.tss.UD500bp.rename.input.bed
depth=1
chr=0
#-----------------------------------------------------------------

#creat step1
if [ ! -d ${path}/shell ]; then
     mkdir -p ${path}/shell
fi
less ${list} |while read a b c ;do echo "sh $shell ${c} $bed $ref $depth $path/${a}/${b}/ ${chr}";done > ${path}/shell/step01_1.get_matrix.sh

#creat step2
echo -e "\
if [ ! -d ${path}/merge ]; then
     mkdir -p ${path}/merge
fi
ls ${path}/*/*/*binDep|head -n 1|while read i;do awk 'BEGIN{FS=OFS=\"\\\t\"}{print \$1,\$2}' \${i};done>${path}/merge/part1.binDep
ls ${path}/*/*/*coverage|head -n 1|while read i;do awk 'BEGIN{FS=OFS=\"\\\t\"}{print \$1,\$2,\$3,\$4}' \${i};done>${path}/merge/part1.coverage
ls ${path}/*/*/*TSS_score|head -n 1|while read i;do awk 'BEGIN{FS=OFS=\"\\\t\"}{print \$1}' \${i};done>${path}/merge/part1.TSS_score
ls ${path}/*/*/*binDep|while read i;do awk '{print \$3}' \${i} >\${i}.1;done
ls ${path}/*/*/*coverage|while read i;do awk '{print \$5}' \${i} >\${i}.1;done
ls ${path}/*/*/*TSS_score|while read i;do awk '{print \$2}' \${i} >\${i}.1;done
echo -e \"gene\\\tbin\$(ls ${path}/*/*/*binDep|while read i;do basename \${i}| awk -F \".\" '{print \$1}';done|while read j;do echo -e \"\\\t\${j}\\\c\";done)\" >${path}/merge/all.binDep
echo -e \"chr\\\tgene\\\tpos\\\tsite\$(ls ${path}/*/*/*coverage|while read i;do basename \${i}| awk -F \".\" '{print \$1}';done|while read j;do echo -e \"\\\t\${j}\\\c\";done)\" >${path}/merge/all.coverage
echo -e \"gene\$(ls ${path}/*/*/*TSS_score|while read i;do basename \${i}| awk -F \".\" '{print \$1}';done|while read j;do echo -e \"\\\t\${j}\\\c\";done)\" >${path}/merge/all.TSS_score
paste ${path}/merge/part1.binDep \$(ls ${path}/*/*/*binDep.1 |while read i;do echo -e \"\\\t\${i}\\\c\";done) >>${path}/merge/all.binDep
paste ${path}/merge/part1.coverage \$(ls ${path}/*/*/*coverage.1 |while read i;do echo -e \"\\\t\${i}\\\c\";done) >>${path}/merge/all.coverage
paste ${path}/merge/part1.TSS_score \$(ls ${path}/*/*/*TSS_score.1 |while read i;do echo -e \"\\\t\${i}\\\c\";done) >>${path}/merge/all.TSS_score
rm -f ${path}/*/*/*binDep.1  ${path}/*/*/*coverage.1 ${path}/*/*/*TSS_score.1 ${path}/merge/part1.binDep ${path}/merge/part1.coverage {path}/merge/part1.TSS_score
">${path}/shell/step02_1.merge_matrix.sh
