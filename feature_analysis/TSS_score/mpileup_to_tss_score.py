# mpileup to depth
import sys
import numpy as np

if len(sys.argv)!=6:
    exit("python3 " + sys.argv[0] + " <bed> <mpileup> <mean_depth> <file_out>\n")


def bed2dic(bed):# 0-based
    with open(bed) as f1:
        bed_dic={}
        for record in f1:
            chr,tssUp,tssDown,gene=record.strip().split('\t')
            bed_dic[gene]=[chr,int(tssUp),int(tssDown)]
        #print(bed_dic)
        return bed_dic


def mplileup2dic(mpileup,mean_depth):
    with open(mpileup) as f2:
        mpileup_dic={}
        for record in f2:
            chr,pos,_,depth,_,_=record.strip().split('\t')
            mpileup_dic[chr+'_'+pos]=int(depth)/float(mean_depth)#######
    #print(mpileup_dic)
    return mpileup_dic


def dic_merge(bed_dic,mpileup_dic,tss_extend):
    gene_dic={}
    for key in bed_dic:
        gene_dic[key]=[None]*(int(tss_extend))
        chr,tssStart,tssEnd=bed_dic[key]
        n=0
        #for pos in range(tssStart+1,tssEnd+1):   #-500,499
        for pos in range(tssStart+1,tssEnd):      #-500,500
            if chr+'_'+str(pos) in mpileup_dic.keys():
                gene_dic[key][n]=mpileup_dic[chr+'_'+str(pos)]
                n=n+1
            else:
                gene_dic[key][n]=0
                n=n+1
    #print(gene_dic)
    return gene_dic

def dic_bin(gene_dic,file_out,tss_extend):
    gene_bin={}
    with open(file_out,'w') as o:
        for key in gene_dic.keys():
            gene_bin[key]=[0]*20
            for num_bin in range(0,20,1):
                gene_bin[key][num_bin]=np.mean(np.array(gene_dic[key][num_bin*50:(num_bin+1)*50])) ########50 should be changed if tss region length is not 1000
            print(key,gene_bin[key],sep='\t',file=o)


def main(bed,mpileup,mean_depth,file_out,tss_extend):
    bed_dic = bed2dic(bed)
    mpileup_dic = mplileup2dic(mpileup,mean_depth)
    gene_dic = dic_merge(bed_dic,mpileup_dic,tss_extend)
    dic_bin(gene_dic,file_out,tss_extend)

main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
#main('/hwfssz5/ST_HEALTH/P20Z10200N0041/lilingguo/PROJECT/HAINAN/TSS/TSS_bed/chr22_500.bed','/hwfssz5/ST_HEALTH/P20Z10200N0041/lilingguo/PROJECT/HAINAN/TSS/result/Control2/depth/Control2_chr22.mplieup',5.1,'/hwfssz5/ST_HEALTH/P20Z10200N0041/lilingguo/PROJECT/HAINAN/TSS/result/Control2_chr22.depth')
