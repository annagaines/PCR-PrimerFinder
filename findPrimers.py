#!/bin/env python3

'''
This is a tool to find primer seqeunces specific to a species or group of speices within  a genus. E.g. A primer specific for Bordetella pertussis
This code takes as input two files, one file contains only the genomes of the targeted sequences. The second files contains that list of all genomes not being targeted.
^wrong

This is a tool to make genus specific primers. E.g. pan-bordetella primers
'''
import getopt 
import sys
import os
import re 
import subprocess
import logging
import shutil
from multiprocessing import Pool
from collections import Counter
from itertools import islice

__group_1_file__ = sys.argv[1]

filelist = [line.rstrip('\n') for line in open(__group_1_file__)]
#print(filelist)

if os.path.isdir("kmers"):
	shutil.rmtree("kmers")
if os.path.isfile("kmers"):
	os.remove("kmers")
os.mkdir("kmers")



def get_kmers(file):
	file = file.rstrip()
	#Jellyfish commands to generate the kmers
	get_kmer_count = f"jellyfish count -C -m 21 -s 100M -o kmers/{(file).replace('_genomic.fna','')}.jf assemblies/{file}"
	get_kmers = f"jellyfish dump kmers/{(file).replace('_genomic.fna','')}.jf |grep -v '>' > kmers/{(file).replace('_genomic.fna','')}.kmers"
	print(get_kmer_count)
	subprocess.check_output(get_kmer_count, shell=True)		
	print(get_kmers)
	subprocess.check_output(get_kmers, shell=True)
	
def combine_freq_and_filt():
	combine_kmers = "cat kmers/*.kmers |sort |uniq -c > combined_uniq_kmers.txt" 
	freq_filt_kmers = " cat combined_uniq_kmers.txt |awk '$1 >= 3 {print  $2}'  |nl |awk '{print \">kmer_\"$1 \"\\n\"$2}' >uniq_kmers_freq_filt.fasta"
	subprocess.check_output(combine_kmers, shell=True)
	subprocess.call(freq_filt_kmers, shell=True)

uniq_kmer_file = open("uniq_kmers_freq_filt.fasta","r")

def calc_tm(file):
	#Read the kmer fasta file in 2 lines at a time, generate dictionary that will save: kmer, kmer seq, qual, Tm, GC
		for line in iter(lambda: list(islice(file,2)),()):
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
					with open("uniq_kmers_freq_bio_filt.fasta", 'w') as out_file:
						out_file.write(f"{line[0]}\n{line[1]}")

#def off_target_check():


def check_dependencies():
	try:
		subprocess.run(["jellyfish"], stdout=devnull, stderr=devnull)
		print("Jellyfish: PASSED")
	except:
		sys.stderr.write("\n WARNING: Cannot find Jellyfish. Check that Jellyfish is downloaded and in your PATH\n\n")
		sys.exit()

	try:
		subprocess.run(['blastn'], stdout=devnull, stderr-devnull)
		print("Blast: PASSED")
	except:
		sys.stderr.write("\n WARNING: Cannot find Jellyfish. Check that Jellyfish is downloaded and in your PATH\n\n")
		sys.exit()



def main():
	with Pool(10) as pool:
		pool.map(get_kmers, filelist)
	combine_and_freq()
	calc_tm(uniq_kmer_file)
	


if __name__ == '__main__':
	main()


