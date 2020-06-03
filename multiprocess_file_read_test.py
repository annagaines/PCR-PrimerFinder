#!/bin/env python3

from multiprocessing import Pool
import sys
import os
from functools import partial
from collections import Counter
from itertools import islice
import subprocess
import re

#g1_uniq_kmers_freq_filt.fasta

def split_kmer_files(file,group):
	split_cmd = f"split -l 500000 {file} {group}_kmers_"
	subprocess.run(split_cmd, shell=True)
	directory = os.getcwd()
	#all_first_reads = [x for x in os.listdir(directory) if "_R1_" in x]
	file_list = [x for x in os.listdir(directory) if "g1_kmers_" in x]
	return file_list

def read_kmer_file(file, group):
	print(f"Filtering kmers from {file}")
	with open(file,'r') as in_file:
		for line in iter(lambda: list(islice(in_file,2)),()):
			line = [l.strip() for l in line]
			if len(line) <2:
				break
			kmer = line[0].replace(">","")
			seq = line[1]
		m = re.search(r'([ACGT])\1{3,}',seq)
		m2 = re.search(r'([CG]){3,}',seq)
		if m or m2:
			pass
		else:
		#Look at sequences to get Tm and GC content	
			l = len(seq)
			counts = Counter(seq.upper())
			#This Tm equation was taken from https://primerdigital.com/fastpcr/m7.html with 2.75 added because it brings the Tm closest to the Tm reported by the tool used by RDB http://www.operon.com/tools/oligo-analysis-tool.aspx
			temp = 69.3 +((41*(counts['G']+counts['C'])-650)/l) 
			if 49 < float(temp) < 66:
				# with open(f"{group}_uniq_kmers_freq_bio_filt.fasta", 'w') as out_file:
				# 	out_file.write(f"{line[0]}\n{line[1]}")
				# out_file.close()
				print(f"{line[0]}\n{line[1]}")

def file_reader(in_file):
	with open(in_file, 'r') as file_in:
		yield islice(file_in,2)

def main():
	# split_kmer_files("g1_uniq_kmers_freq_filt.fasta", 'g1')
	pool = Pool(10)
	infile_list = split_kmer_files("g1_uniq_kmers_freq_filt.fasta", 'g1')
	pool.map(partial(read_kmer_file, group='g1'),infile_list)			
	


if __name__ == '__main__':
	main()
