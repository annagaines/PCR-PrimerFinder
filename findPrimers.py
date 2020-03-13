#!/bin/env python3

'''
This is a tool to find primer seqeunces specific to a species or group of speices within  a genus. E.g. A primer specific for Bordetella pertussis
This code takes as input two files, one file contains only the genomes of the targeted sequences. The second files contains that list of all genomes not being targeted.

'''
import getopt 
import sys
import os
import re 
import subprocess
import logging
import shutil
from multiprocessing import Pool

__group_1_file__ = sys.argv[1]

filelist = [line.rstrip('\n') for line in open(__group_1_file__)]
print(filelist)

def get_kmers(file):
	file = file.rstrip()
	get_kmer_count = f"jellyfish count -C -m 21 -s 100M -o kmers/{(file).replace('_genomic.fna','')}.jf assemblies/{file}"
	get_kmers = f"jellyfish dump -s 100M -o kmers/{(file).replace('_genomic.fna','')}.jf assemblies/{file}"
	print(get_kmer_count)
	subprocess.check_output(get_kmer_cmd, shell=True)		

# def get_genome_kmers(__group_1_file__):
# 	if os.path.isdir("kmers"):
# 		shutil.rmtree("kmers")	
# 	os.mkdir(("kmers"))
# 	for line in __group_1_file__:
	


def main():
	with Pool(10) as pool:
		pool.map(get_kmers, filelist)


	#get_genome_kmers(__group_1_file__)

if __name__ == '__main__':
	main()


