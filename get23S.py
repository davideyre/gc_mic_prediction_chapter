#!/usr/bin/env python
import os, pysam, sys

def get23S(bwaPath, inputStem, outputStem, geneLibPath, sampleName, outFile):

	#bwaPath - path to bwa mem
	#inputStem - prefix for input fastq files which are assumed to end R1.fastq and R2.fastq for forward and reverse reads
	#outputStem - prefix for intermediate output file
	#geneLibPath - location of reference_23s.fa
	#sampleName - human readable sample name
	#outfile - name for output file
	
	o = open(outFile, 'w')
	header = 'sample_name\tantibiotic\tgene\tsite\taa_snp_presence\n'
	o.write(header)
	
	ref = '%s/reference_23s.fa'%geneLibPath
	fq1 = inputStem+'.R1.fastq'
	fq2 = inputStem+'.R2.fastq'
	
	out_sam = inputStem+'.sam'
	out_bam = inputStem+'.bam'
	
	#map files using bwa mem
	cmd = '%sbwa mem -t 4 %s %s %s > %s'%(bwaPath, ref, fq1, fq2, out_sam)
	print cmd
	os.system(cmd)
	
	# put the mapped reads (-F 4) into sorted Bam
	cmd = 'samtools view -uS -F 4 %s | samtools sort - -o %s'%(out_sam, out_bam)
	print cmd
	os.system(cmd)
	
	cmd = 'samtools index %s'%out_bam
	print cmd
	os.system(cmd)
	
	#get A2059G == 2045 numbered starting from 1
	#get C2611T == 2597 numbered starting from 1
	base2045 = []
	base2597 = []
	
	samfile = pysam.Samfile(out_bam, "rb")
	for pileupcolumn in samfile.pileup( '23s_rna', 2000, 2600):
		for pileupread in pileupcolumn.pileups:
			if pileupcolumn.pos == 2044 and pileupread.query_position != None:
				base = pileupread.alignment.query_sequence[pileupread.query_position]
				quality = ord(pileupread.alignment.qual[pileupread.query_position])
				if quality-33>=30:
					base2045.append(base)
			elif pileupcolumn.pos == 2596 and pileupread.query_position != None:
				base = pileupread.alignment.query_sequence[pileupread.query_position]
				quality = ord(pileupread.alignment.qual[pileupread.query_position])
				if quality-33>=30:
					base2597.append(base)
	
	outputLog = outputStem+'_azt_bwa_log.txt'

	ol = open(outputLog, 'w')
	header = 'guid\tsite\tA\tC\tG\tT\n'
	ol.write(header)
	
	a = len([b for b in base2045 if b=='A'])
	c = len([b for b in base2045 if b=='C'])
	g = len([b for b in base2045 if b=='G'])
	t = len([b for b in base2045 if b=='T'])
	ol.write('%s\t2045\t%s\t%s\t%s\t%s\n'%(sampleName, a, c, g, t))
	
	a = len([b for b in base2597 if b=='A'])
	c = len([b for b in base2597 if b=='C'])
	g = len([b for b in base2597 if b=='G'])
	t = len([b for b in base2597 if b=='T'])
	ol.write('%s\t2597\t%s\t%s\t%s\t%s\n'%(sampleName, a, c, g, t))
	
	ol.close()
	
	#read in copy number for A2059G and C2611T, aka 2045 and 2597
	f = open(outputLog, 'r')
	header = f.next() #guid	site	A	C	G	T
	baseCountDict = dict()
	
	for l in f:
		l = l.strip().split()
		base = int(l[1])
		if base == 2045:
			snp_column = 4
		elif base == 2597:
			snp_column = 5
		propSnp = float(int(l[snp_column])) / sum(int(i) for i in l[2:]) #check what proportion T - the SNP for both
		copyNumber = int(round(propSnp * 4))
		baseCountDict[base] = copyNumber	
	
	f.close()
	
	for base in baseCountDict.keys():
		if base==2045:
			mut = 'A2059G'
		else:
			mut = 'C2611T'
		outputString = '%s\t%s\t%s\t%s\t%s\n'%(sampleName, 'azt', 'rrna_mutation_copy_number', mut, baseCountDict[base])
		sys.stdout.write(outputString)
		o.write(outputString)