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

__group_1_file__ = open(sys.argv[1],"r")
__group_2_file__ = open(sys.argv[2], "r")

group1 = sys.argv[1].replace(".txt",'')
group2 = sys.argv[2].replace(".txt",'')

def get_genome_kmers(__group_1_file__,__group_2_file__, group1, group2):
	if os.path.isdir(group1):
		shutil.rmtree(group1)	
	os.mkdir((group1).replace('.txt',''))
	for line in __group_1_file__:
		line = line.rstrip()
		get_kmer_count = f"jellyfish count -C -m 21 -s 100M -o {group1}/{(line).replace('_genomic.fna','')}.jf assemblies/{line}"
		get_kmers = f"jellyfish dump -C -m 21 -s 100M -o {group1}/{(line).replace('_genomic.fna','')}.jf assemblies/{line}"
		#subprocess.check_output(get_kmer_cmd, shell=True)		
	if os.path.isdir(group2):
		shutil.rmtree(group2)	
	os.mkdir((group2).replace('.txt',''))
	for line in __group_2_file__:
		line = line.rstrip()
		get_kmer_count = f"jellyfish count -C -m 21 -s 100M -o {group2}/{(line).replace('_genomic.fna','')}.jf assemblies/{line}"
		get_kmers = f"jellyfish dump -C -m 21 -s 100M -o {group2}/{(line).replace('_genomic.fna','')}.jf assemblies/{line}"
		#subprocess.check_output(get_kmer_cmd, shell=True)

# def kmerize_genomes():


# def split_genomes():


def main():
	get_genome_kmers(__group_1_file__,__group_2_file__, group1, group2)

if __name__ == '__main__':
	main()


