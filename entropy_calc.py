#!/bin/env python3

import sys
import re
from math import log
from collections import Counter

#This script is used to calculates the entopy of the provided kmer and matches them with their good or bad designation
#./entropy_calc.py round1.diffKmers.fa blastValPass.txt

file = open(sys.argv[1],"r")
file2 = sys.argv[2]

good_kmers = [line.rstrip('\n') for line in open(file2)]


def calc_entropy(file):
	for line in file:
		line = line.rstrip()
		if line.startswith('>'):
			kmer = line.replace(">","")
			if kmer in good_kmers:
				print(kmer,"good")
			else:
				print(kmer, "bad")
		else:
			l = len(line)
			counts = Counter(line)
			parts = [(i/l)*(log((i/l),2)) if i != 0 else 0 for i in counts.values()]
			Shannon = -sum(parts)
			print(Shannon)

def main():
	calc_entropy(file)

if __name__ == '__main__':
	main()

