import sys
import pandas as pd
import glob
import operator
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


"""
#NC_045512.2:21563-25384	1	A	1	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:1:42.00:17.00:0.00:1:0:0.00:0.02:0.00:0:0.00:0.00:0.
00	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0
.00:0.00	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
"""




import sys
import pandas as pd
import glob
import operator


files=glob.glob("*.nucl_out")
base={}

for f in files:
	samp=f.strip().split("/")[-1].split(".")[0]
	inF=open(f.strip(),"r")
	
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
				

				if LETTER==M and PERC>=0.01: #MAGGIORITARI
						
					if samp not in base.keys():
						base[samp]={}
					base[samp][pos]=LETTER
						
hash={}
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
				

				if LETTER!=M and PERC>=0.01: # DO STUFF!!!!!!!
					codpos=((pos-301)%3)+1
					if codpos==1:
						CODorig=base[samp][pos]+base[samp][pos+1]+base[samp][pos+2]
						CODmut=LETTER+base[samp][pos+1]+base[samp][pos+2]
					if codpos==2:
						CODorig=base[samp][pos-1]+base[samp][pos]+base[samp][pos+1]
						CODmut=base[samp][pos-1]+LETTER+base[samp][pos+1]
					if codpos==3:
						CODorig=base[samp][pos-2]+base[samp][pos-1]+base[samp][pos]
						CODmut=base[samp][pos-2]+base[samp][pos-1]+LETTER
					MUTtype="NA"
					if "In" not in CODorig and "Del" not in CODorig:
						if "In" in CODmut or "Del" in CODmut:
							MUTtype="Indel"
							AAorig="NA"
							AAmut="NA"
						else:
							AAorig=Seq(CODorig,generic_dna).translate()
							AAmut=Seq(CODmut,generic_dna).translate()
							if AAorig == AAmut:
								MUTtype="Syn"
							elif AAmut=="*":
								MUTtype="STOP"
							else:
								MUTtype="NotSyn"

					#print CODorig, CODmut, AAorig, AAmut, MUTtype
					if MUTtype not in hash[samp].keys():
						hash[samp][MUTtype]=0
					hash[samp][MUTtype]+=1





TAB=pd.DataFrame.from_dict(hash, orient="index")
#print TAB

TAB.to_csv("CountSynNotSyn.csv",sep=";")
				

