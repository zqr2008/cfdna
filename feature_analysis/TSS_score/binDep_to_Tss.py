import numpy as np
import sys

fi=sys.argv[1]
fo=sys.argv[2]

with open(fi,'r') as f, open(fo,'w') as o:
    #for line in f:
    #    column_names=[]
    #    column_names.append(line.strip().split('\t')[0])
    #    column_names=column_names+line.strip().split('\t')[2:]
    #    col="\t".join(column_names)
    #    print(col,file=o)
    #    break
    n=0
    for line in f:
        gene=line.strip().split('\t')[0]
        if n==0 :
            if 'None' in line:
                print(gene,'',sep='\t',file=o)
                n=n+1
            else:
                values=list(map(float,line.strip().split('\t')[2:]))
                matrix=np.array(values)
                n=n+1
        elif n<9 :
            if 'None' in line:
                n=n+1
            else:
                values=list(map(float,line.strip().split('\t')[2:]))
                matrix=matrix+np.array(values)
                n=n+1
        elif n==9:
            if 'None' in line:
                n=0
            else:
                values=list(map(float,line.strip().split('\t')[2:]))
                matrix=(matrix+np.array(values))/10
                col="\t".join(list(map(str,list(matrix))))
                print(gene,col,sep='\t',file=o)
                n=0
