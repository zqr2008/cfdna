#mpileup to per site coverage

import sys

def mplileup2dic(mpileup,mean_depth):
    with open(mpileup) as f2:
        mpileup_dic={}
        for record in f2:
            chr,pos,_,depth,_,_=record.strip().split('\t')
            mpileup_dic[chr+'_'+pos]=int(depth)/float(mean_depth)
    return mpileup_dic


def bedProcess(line):
            chr,TssStart,TssEnd,gene=line.strip().split('\t')
            return chr,int(TssStart),int(TssEnd),gene

def main_process(chr,TssStart,TssEnd,gene,mpileup_dic,o):
        mylist=range(TssStart+1,TssEnd+1,1)
        site=-500
        for pos in mylist:
            if chr+'_'+str(pos) in mpileup_dic.keys():
                print(chr,gene,pos,site,mpileup_dic[chr+'_'+str(pos)],sep='\t',file=o)
                site=site+1
            else:
                print(chr,gene,pos,site,'0',sep='\t',file=o)
                site=site+1

def main(bed,mpileup,mean_depth,file_out):
    with open(bed) as f1 , open(file_out,'w') as o:
        #print('chr','gene','pos','site','coverage',sep='\t',file=o)
        mpileup_dic=mplileup2dic(mpileup,mean_depth)
        for line in f1:
            chr,TssStart,TssEnd,gene=bedProcess(line)
            main_process(chr,TssStart,TssEnd,gene,mpileup_dic,o)

main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

