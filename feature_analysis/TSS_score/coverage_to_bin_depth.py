import sys
import numpy as np

def main(coveragefile,binfile):
    with open(coveragefile) as fi, open(binfile,'w') as fo:
        dic={}
        for line in fi:
            chr,gene,pos,site,coverage=line.strip().split('\t')
            if gene not in dic.keys():
                dic[gene]=[float(coverage)]
            else:
                dic[gene].append(float(coverage))
        for key in dic.keys():
            n=0
            depList=dic[key]
            flanking=np.mean(np.array(depList[0:250]+depList[750:1000]))
            for bin in range(0,10):
                start=250+bin*50
                end=300+bin*50
                if flanking!=0:
                    binCov=np.mean(np.array(depList[start:end]))/flanking
                    print(key,'bin'+str(bin),binCov,sep='\t',file=fo)
                else:
                    print(key,'bin'+str(bin),None,sep='\t',file=fo)

main(sys.argv[1],sys.argv[2])
