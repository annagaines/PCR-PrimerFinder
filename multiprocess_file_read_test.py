#!/bin/env python3

from multiprocessing import Pool
import sys
from functools import partial
from collections import Counter
from itertools import islice
import re

#g1_uniq_kmers_freq_filt.fasta

def read_kmer_file(lines, group):
		print(lines.rstrip())
		if lines.startswith('>'):	
			kmer = lines[0].replace(">","") 
		else:
			seq = lines[1]
		print(kmer)
		print(seq)

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
	pool = Pool(10)
	infile = "g1_uniq_kmers_freq_filt.fasta"
	lines = file_reader(infile)
	pool.map(partial(read_kmer_file, group='g1'),next(lines))			
	#infile.close()


if __name__ == '__main__':
	main()
