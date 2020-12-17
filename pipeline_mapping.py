import sys
import os

#os.system("bowtie2-build Reference_Spike+300.fasta REFERENCE")
#os.system("rm -rf Depthplots")
#os.system("mkdir Depthplots")

inF=open(sys.argv[1].strip(),"r")

for i in inF.xreadlines():
	i=i.strip().split(",")
	FILE=i[0]
	NAME=i[1]
	cmd="fastqc %s1_001.fastq.gz %s2_001.fastq.gz" %(FILE,FILE)
	os.system(cmd)
	cmd="./fastp --in1 %s1_001.fastq.gz --in2 %s2_001.fastq.gz --out1 %s_trim1.fastq.gz --out2 %s_trim2.fastq.gz -f 26 -F 26 -t 28 -T 28 -l 100 --average_qual 25" %(FILE,FILE,NAME,NAME)
	os.system(cmd)
	cmd="bowtie2 --very-sensitive -p 11 -x REFERENCE -1 %s_trim1.fastq.gz -2 %s_trim2.fastq.gz -S %s.sam" %(NAME,NAME,NAME)
	os.system(cmd)
	cmd="echo %s" %(NAME)
	os.system(cmd)
	cmd='grep -v "^@" %s.sam | cut -f 3 | sort | uniq -c' %(NAME)
	os.system(cmd)
	cmd="/usr/bin/samtools view -S -b %s.sam > %s.bam" %(NAME,NAME)
	os.system(cmd)
    cmd="/usr/bin/samtools sort %s.bam %s.sorted" %(NAME,NAME)
    #os.system(cmd)
    #cmd="/usr/bin/samtools mpileup -f Reference_Spike+300.fasta %s.sorted.bam | cut -f2,4 | awk '$1>300' | awk '$1<4123' | awk '$1=$1-300' > %s.depth_each_base" %(NAME,NAME)
    #os.system(cmd)
	#cmd="Rscript Plot_coverage.R %s.depth_each_base > %s_depthlog.txt" %(NAME,NAME)
	#os.system(cmd)
	#os.system("mv *.png Depthplots/")
    cmd="./bam-readcount -w 0 -i -f Reference_Spike+300.fasta %s.sorted.bam > %s.nucl_out" %(NAME,NAME)
	os.system(cmd)
	cmd="java -XX:+UseParallelGC -XX:NewRatio=9 -Xms7G -Xmx20G -jar clique-snv.jar -m snv-illumina -in %s.sam -tf 0.01 -t 3 -cm 'fast' -outDir ./%s_clique/" %(NAME,NAME)
	os.system(cmd)
	cmd="python trim_haplotypes.py %s_clique/%s.fasta Reference_spike_gene.fasta" %(NAME,NAME)
	os.system(cmd)

inF.close()

os.system("python format_readcount2table.py")
os.system("python format_mutation_patterns.py")
os.system("python mutation_rate2.py")
os.system("python fasta_haplotype_count.py")
os.system("python 1_STOP_FRAME_count.py")


