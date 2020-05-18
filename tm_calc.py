#!/bin/env python3


#!/bin/env python3

import sys
import re
from collections import Counter
from itertools import islice

#This script is used to calculate the melting temperature of the kmers nad matches them with their good or bad designation using blastValPass.txt (which are the kmers with High IMportance in synth.valid.tsv)
#./tm_calc.py round1.diffKmers.fa blastValPass.txt

file = open(sys.argv[1],"r")
file2 = sys.argv[2]

good_kmers = [line.rstrip('\n') for line in open(file2)]

kmers = {}

def calc_tm(file):
	#Read the kmer fasta file in 2 lines at a time, generate dictionary that will save: kmer, kmer seq, qual, Tm, GC
	for line in iter(lambda: list(islice(file,2)),()):
		line = [l.strip() for l in line]
		if len(line) <2:
			break
		kmer = line[0].replace(">","")
		seq = line[1]
		if kmer not in kmers:
			kmers[kmer] = {'seq' : seq, 'qual' : '','gc' : 0, 'tm' : 0 }
	#Compare kmers to list of good kmers, assign quality
	for kmer in kmers:
		if kmer in good_kmers:
			kmers[kmer]['qual'] = 'good'
		else:
			kmers[kmer]['qual'] = 'bad'
	#Look at sequences to get Tm and GC content
		l = len(kmers[kmer]['seq'])
		counts = Counter(kmers[kmer]['seq'].upper())
		#This Tm equation was taken from https://primerdigital.com/fastpcr/m7.html with 2.75 added because it brings the Tm closest to the Tm reported by the tool used by RDB http://www.operon.com/tools/oligo-analysis-tool.aspx
		temp1 = 69.3 +((41*(counts['G']+counts['C'])-650)/l) 
		temp2 = (4*(counts['G']+counts['C'])) + (2*(counts['A']+counts['T']))
		print("-------------------")
		print(f"Temp1: {temp1}")
		print(f"Temp2: {temp2}")

		kmers[kmer]['tm'] = "{:.2f}".format(temp1)
		gc = (counts['G']+counts['C'])*100/(counts['G']+counts['C']+counts['A']+counts['T'])
		kmers[kmer]['gc'] = "{:.2f}".format(gc)	
		
def summarize(kmer_dict):
	pass_list = []
	gc_pass_count = 0
	tm_pass_count = 0
	for kmer in kmer_dict:
		# if 40 < float(kmer_dict[kmer]['gc']) < 60:
		# 	gc_pass_count += 1
		if 57 < float(kmer_dict[kmer]['tm']) < 64:
			tm_pass_count += 1
			pass_list.append(kmer)
	print(f"Number passed: {tm_pass_count}")
	print(f"Number passed: {tm_pass_count*100/len(kmer_dict)}")

def main():
	calc_tm(file)
	summarize(kmers)
if __name__ == '__main__':
	main()

