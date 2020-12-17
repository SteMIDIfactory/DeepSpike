
#NC_045512.2:21563-25384	1	A	1	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:1:42.00:17.00:0.00:1:0:0.00:0.02:0.00:0:0.00:0.00:0.00	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00

import sys
import pandas as pd
import glob
import operator

files=glob.glob("*.nucl_out")
hash={}
for f in files:
	samp=f.strip().split("/")[-1].split(".")[0]
	inF=open(f.strip(),"r")

	hash[samp]={}
	for i in inF.xreadlines():
		i=i.strip().split("\t")
		pos=int(i[1])
		#print pos
		if pos>=301 and pos <=4122:
			tot=0
			for x in range(4,len(i)):
				tot=tot+int(i[x].strip().split(":")[1])
			CONS={}
			for x in range(4,len(i)):
				LETTER=i[x].strip().split(":")[0][0]
				if LETTER=="+":
					LETTER="In"
				if LETTER=="-":
					LETTER="Del"
				VALUE=i[x].strip().split(":")[1]
				PERC=float(VALUE)/float(tot)
				CONS[LETTER]=PERC
			M=max(CONS.iteritems(), key=operator.itemgetter(1))[0]
			for x in range(4,len(i)):
				LETTER=i[x].strip().split(":")[0][0]
				if LETTER=="+":
					LETTER="In"
				if LETTER=="-":
					LETTER="Del"
				VALUE=i[x].strip().split(":")[1]
				PERC=float(VALUE)/float(tot)
				CONS[LETTER]=PERC
				if LETTER!=M and PERC>=0.01:
					MUT="%s_%s" %(M,LETTER)
					if MUT not in hash[samp].keys():
						hash[samp][MUT]=0
					hash[samp][MUT]+=1

	inF.close()



TAB=pd.DataFrame.from_dict(hash, orient="columns")

TAB.to_csv("MutationPatterns.csv",sep=";")
