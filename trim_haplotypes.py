import sys
from Bio import SeqIO
import os
import subprocess

hap=sys.argv[1].strip()
inF=open(hap,"r")

count=0
for i in SeqIO.parse(inF,"fasta"):
	count+=1

#print count


ref=sys.argv[2].strip()

cmd="makeblastdb -in %s -dbtype nucl -logfile log" %(hap)
os.system(cmd)


cmd="blastn -query %s -db %s -outfmt '6 sseqid qlen length pident sseq' | awk '$4>90' | awk '$3<$2+10' | awk '$3>$2-10'" %(ref,hap)
output = subprocess.check_output(cmd, shell=True)
output=output.strip().split("\n")

hash={}
for x in output:
	x=x.split()
	prevalence=float(x[0].strip().split("_")[2])
	if x[4] not in hash.keys():
		hash[x[4]]=0.0
	hash[x[4]]+=prevalence

import operator
MM=max(hash.iteritems(), key=operator.itemgetter(1))[0]
#print MM


outname=hap.strip().split("/")[-1].split(".")[0]+"_haplotypes.fasta"

recordS=[]
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
name="%s___BASE___%s" %(hap.strip().split("/")[-1].split(".")[0], str(round(hash[MM],3)))
record = SeqRecord(Seq(MM,IUPAC.ambiguous_dna),id=name, description="")
recordS.append(record)

#SeqIO.write(record,outname,"fasta")

del hash[MM]

SO=sorted(hash.iteritems(), key=operator.itemgetter(1))
SO.reverse()
#print SO
c=0
for y in SO:
	c+=1
	name="%s___%i___%s" %(hap.strip().split("/")[-1].split(".")[0],c, str(round(y[1],3)))
	record = SeqRecord(Seq(y[0],IUPAC.ambiguous_dna),id=name, description="")
	recordS.append(record)
SeqIO.write(recordS,outname,"fasta")

#print record


inF.close()
