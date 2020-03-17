#!/bin/env python3

import sys
import re
from math import log
from collections import Counter

#This script is used to calculate the accuracy, precision, specificity, and sensitivity of entropy cut off

file = open(sys.argv[1],"r")

kmer_entropy = {}

def calc_entropy(file):
	for line in file:
		a_count = 0
		t_count = 0
		g_count = 0
		c_count = 0
		line = line.rstrip()
		if line.startswith('>'):
			kmer = line.replace(">","")
			if kmer not in kmer_entropy:
				kmer_entropy[kmer] = {}
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

