import sys
import pandas as pd
import glob
import operator

regions={"NTD":[[37,915]],"RBD":[[955,1623]],"SD1":[[1624,2055]],"SD2":[[2056,2361]],"FP":[[2362,2418]],"HR1":[[2734,2952]],"HR2":[[3406,3639]],"TM":[[3640,3711]],"S1":[[1,2055]],"S2":[[2056,3822]]}

import sys
import pandas as pd
import glob
import operator


files=glob.glob("*.nucl_out")
hash={}
counter={}
for f in files:
	samp=f.strip().split("/")[-1].split(".")[0]
	inF=open(f.strip(),"r")
	
	hash[samp]={}
	for i in inF.xreadlines():
		i=i.strip().split("\t")
		pos=int(i[1])
		
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
				#creo hash in cui il counter e zero
				
				if LETTER!=M and PERC>=0.01: #se if e true il conteggio aumenta di 1 perche la mut nn era nel base
					MUT="%s_%s" %(M,LETTER)
					if samp not in counter.keys():				
						counter[samp]=0
					counter[samp]+=1
					for w in regions.keys():
						for z in regions[w]:
							if pos in range(z[0],z[1]+1):
								if w not in hash[samp].keys():
									hash[samp][w]=0
								hash[samp][w]+=1
											
					

	inF.close()
	



TAB=pd.DataFrame.from_dict(hash, orient="index")

TAB.to_csv("MutationCount_Regions_new.csv",sep=";")				


