
#NC_045512.2:21563-25384	1	A	1	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:1:42.00:17.00:0.00:1:0:0.00:0.02:0.00:0:0.00:0.00:0.00	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00

import sys
import pandas as pd
import glob
import operator
from Bio import SeqIO

files=glob.glob("*.nucl_out")

hash={}
consensusDICT={}
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
			for y in ["A","T","C","G","N","In","Del"]:
				ind="%i_%s" %(pos,y)
				hash[samp][ind]=0
			#print hash
			#break
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
				if LETTER!="=" and PERC>=0.01:
        	                	ind="%i_%s" %(pos,LETTER)
					hash[samp][ind]=1
			M=max(CONS.iteritems(), key=operator.itemgetter(1))[0]
			if pos not in consensusDICT.keys():
				consensusDICT[pos]={}
			if M not in consensusDICT[pos].keys():
				consensusDICT[pos][M]=0
			consensusDICT[pos][M]+=1
	inF.close()

#print consensusDICT

consensus={}
for j in consensusDICT.keys():
	best=max(consensusDICT[j].iteritems(), key=operator.itemgetter(1))[0]
	consensus[j]=best



print ">consensus"
c=""
for j in consensus.keys():
	c=c+consensus[j]
print c

refF=open(sys.argv[1].strip(),"r")

for i in SeqIO.parse(refF,"fasta"):
	SEQ=str(i.seq).upper()
	print ">Wuhan"
	print SEQ



refF.close()

Wuhan={}
for x in range(0,len(SEQ)):
	Wuhan[x+301]=SEQ[x]


hashCODONmaj={}
hashmaj={}


for f in files:
        samp=f.strip().split("/")[-1].split(".")[0]
        inF=open(f.strip(),"r")
	hashCODONmaj[samp]={}
        hashmaj[samp]={}
	cc=""
        for i in inF.xreadlines():
                i=i.strip().split("\t")
                pos=int(i[1])
                #print pos
                if pos>=301 and pos <=4122:
                        tot=0
                        for x in range(4,len(i)):
                                tot=tot+int(i[x].strip().split(":")[1])
                        #print hash
                        #break
			codon=((pos-301)/3)+1
			if codon not in hashCODONmaj[samp].keys():
				hashCODONmaj[samp][codon]=0
			hashmaj[samp][pos]=0
			LETTERS={}
                        for x in range(4,len(i)):
                                LETTER=i[x].strip().split(":")[0][0]
                                if LETTER=="+":
                                        LETTER="In"
                                if LETTER=="-":
                                       LETTER="Del"
                                VALUE=i[x].strip().split(":")[1]
                                PERC=float(VALUE)/float(tot)
                                CONS[LETTER]=PERC
                                if LETTER!="=" and PERC>=0.01:
					LETTERS[LETTER]=PERC
			M=max(LETTERS.iteritems(), key=operator.itemgetter(1))[0]
			cc=cc+M
			if M!=Wuhan[pos]:
				hashmaj[samp][pos]=1
				hashCODONmaj[samp][codon]=1

	print ">"+samp
	print cc
	inF.close()

TAB=pd.DataFrame.from_dict(hashmaj, orient="columns")

TAB.to_csv("MutatedPosition_maj_wuhan.csv",sep=";")


TAB=pd.DataFrame.from_dict(hashCODONmaj, orient="columns")

TAB.to_csv("MutatedCodon_maj_wuhan.csv",sep=";")





