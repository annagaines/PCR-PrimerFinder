#!/bin/env python3

'''
This is a tool to find primer seqeunces specific to a species or group of speices within  a genus. E.g. A primer specific for Bordetella pertussis
This code takes as input two files, one file contains only the genomes of the targeted sequences. The second files contains that list of all genomes not being targeted.
^wrong

This is a tool to make genus specific primers. E.g. pan-bordetella primers
'''
import click
import sys
import os
import re 
import subprocess
import logging
import shutil
from multiprocessing import Pool
from collections import Counter
from itertools import islice
from functools import partial

uniq_kmer_file = ""
uniq_filtered_kmers = ""
nt_blast_results = ""

@click.command()
#@click.option('-m','--mode','mode', help="Desired mode for primer prediction.'p' for pan-genus prediction (provide assemblies for genomes in this genus and file list of all genomes). 's' for species prediction (provide assemblies for genomes in this genus and two filelist, one for the speices and one with genus other than target species. 'd' for differentiating between two species. (Provide assemblies for both species, supply a file list for each species)")
@click.option('-1', '--input1','genomelist1', help="File with list of target genomes. This is the filename only and does not include the path to where the files are located.All genomes should be located in current working directory in assemblies/" , type=click.File(mode='r'))
@click.option('-2', '--input2','genomelist2', help="File with list of target genomes. This is the filename only and does not include the path to where the files are located. All genomes should be located in current working directory in assemblies/" , type=click.File(mode='r'))
@click.option('-w','--whiteList','whiteList', help="File containing taxids categorized as 'Good' and 'Acceptable' target matches for the identified k-mers",type=click.File(mode='r'))
@click.option('-b','--blackList','blackList', help="File containing taxids categorized as 'Bad' target matches for the identified k-mers",type=click.File(mode='r'))
@click.option('--log', 'log', default=None, type=click.File(mode='w'))
@click.option('-pident','percentIdentity',default=100, help="The minimum percent identity for k-mer match against NCBI nucleotide database")
@click.option('-f','--frequency','frequency' , default=.9 ,help='Desired frequency of k-mer in target genomes. (float value)')
def main(genomelist1,genomelist2, whiteList,blackList,percentIdentity, log, frequency):
	#Check is a kmer directory exists, if it does delete it, if it doesn't make one
	if os.path.isdir("tmp"):
		shutil.rmtree("tmp")
	if os.path.isfile("tmp"):
		os.remove("tmp")
	os.mkdir("tmp")

	if os.path.isdir("kmers"):
		shutil.rmtree("kmers")
	if os.path.isfile("kmers"):
		os.remove("kmers")
	os.mkdir("kmers")

	#create log file
	if log is None: log = "primerFinder.log"
	logging.basicConfig(filename=log, filemode='w', level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	sys.stderr.write(f"\nWriting log file to: {log}\n")
	logging.debug(f"Command: {' '.join(sys.argv)}")
	#if mode == 's':
	#read in the filelist of target genomes
	filelist1 = [line.rstrip('\n') for line in genomelist1]
	filelist2 = [line.rstrip('\n') for line in genomelist2]
	kmer_frequency_1 = len(filelist1) * frequency 
	kmer_frequency_2 = len(filelist2) * frequency 

	group1_genomes, group2_genomes = get_genome_groups(filelist1, filelist2)
	make_jellyfish_directories("g1")
	make_jellyfish_directories("g2")
	with Pool(10) as pool:
		pool.map(partial(get_kmers, group='g1'), filelist1)
		pool.map(partial(get_kmers, group='g2'), filelist2)
	combine_freq_and_filt(kmer_frequency_1, 'g1')
	combine_freq_and_filt(kmer_frequency_2, 'g2')
	file_list1 = split_kmer_files("tmp/g1_uniq_kmers_freq_filt.fasta", 'g1')
	file_list2 = split_kmer_files("tmp/g2_uniq_kmers_freq_filt.fasta", 'g2')
	with Pool(10) as pool:	
		pool.map(partial(calc_tm, group='g1'), file_list1)
		pool.map(partial(calc_tm, group='g2'), file_list2)
	
	compare_kmers("g1","g2")
	blast_against_genomes()

	#off_target_check("g1_uniq_kmers_freq_bio_filt.fasta",percentIdentity)
	# check_nt_blast_results(whiteList,blackList,nt_blast_results)

def get_genome_groups(filelist1, filelist2):
	logging.debug("\tCreating lists of genomes for g1 and g2.")
	group1_genomes = {}
	group2_genomes = {}
	for file in filelist1:
		get_genome_cmd = f"head -n 1 assemblies/{file}|cut -d' ' -f1 |sed 's/>//' "
		genome = subprocess.check_output(get_genome_cmd, shell=True, universal_newlines=True)
		group1_genomes[genome] = "g1"

	for file in filelist2:
		get_genome_cmd = f"head -n 1 assemblies/{file}|cut -d' ' -f1 |sed 's/>//' "
		genome = subprocess.check_output(get_genome_cmd, shell=True, universal_newlines=True)
		group2_genomes[genome] = "g2"
	return	(group1_genomes, group2_genomes)

def make_jellyfish_directories(group):
	logging.debug(f"\tCreating directory for k-mer in {group}: kmers/{group}")
	if os.path.isdir(f"kmers/{group}"):
		shutil.rmtree(f"kmers/{group}")
	if os.path.isfile(f"kmers/{group}"):
		os.remove(f"kmers/{group}")
	os.mkdir(f"kmers/{group}")

#Take each file in filelist and break create .jf and then .kmer files
def get_kmers(file, group):
	logging.debug(f"\tCreating k-mers for {group}: {file}")
	file = file.rstrip()	
	#Jellyfish commands to generate the kmers
	get_kmer_count = f"jellyfish count -C -m 21 -s 100M -o kmers/{group}/{(file).replace('_genomic.fna','')}.jf assemblies/{file}"
	get_kmers = f"jellyfish dump kmers/{group}/{(file).replace('_genomic.fna','')}.jf |grep -v '>' > kmers/{group}/{(file).replace('_genomic.fna','')}.kmers"
	#print(get_kmer_count)
	subprocess.check_output(get_kmer_count, shell=True)		
	#print(get_kmers)
	subprocess.check_output(get_kmers, shell=True)

#Combine seperate .kmer files and filter it based on the given frequency	
def combine_freq_and_filt(frequency, group):
	logging.debug(f"\tFiltering k-mers in {group} using a frequency threshold of {frequency}")
	combine_kmers = f"cat kmers/{group}/*.kmers |sort |uniq -c > tmp/{group}_combined_uniq_kmers.txt" 
	freq_filt_kmers = f"cat tmp/{group}_combined_uniq_kmers.txt |awk '$1 >= {frequency} {{print  $2}}'  |nl |awk '{{print \">kmer_\"$1 \"\\n\"$2}}' > tmp/{group}_uniq_kmers_freq_filt.fasta"
	subprocess.check_output(combine_kmers, shell=True)
	subprocess.call(freq_filt_kmers, shell=True)

def split_kmer_files(file,group):
	logging.debug(f"Spliting kmers from {file} in to files prefixed with {group}_kmers_")
	split_cmd = f"split -l 500000 {file} tmp/{group}_kmers_"
	subprocess.run(split_cmd, shell=True)
	#directory = os.getcwd()
	file_list = [x for x in os.listdir("tmp") if "g1_kmers_" in x]
	return file_list

def calc_tm(file, group):
	#Read the kmer fasta file in 2 lines at a time, generate dictionary that will save: kmer, kmer seq, qual, Tm, GC
	logging.debug(f"\tFiltering k-mers in tmp/{file} using a Tm range of 50 to 66 and removing k-mers with homopolymers > 3 or repeats of G and/or C > 3. ")
	
	with open(f"tmp/{file}",'r') as in_file:
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
					with open(f"tmp/{group}_uniq_kmers_freq_bio_filt.fasta", 'w') as out_file:
						out_file.write(f"{line[0]}\n{line[1]}")
					out_file.close()


def compare_kmers(group1, group2):
	logging.debug(f"\tCreating list of uniq different k-mers.")
	combine_kmer_cmd = f"cat tmp/{group1}_uniq_kmers_freq_bio_filt.fasta tmp/{group2}_uniq_kmers_freq_bio_filt.fasta| sort --parallel 20 -S 10G| uniq -c | awk '$1 == 1 {{print $2}}' |nl | awk '{{print \">kmer_\"$1\"\\n\" $2 }}' > tmp/diffKmers.fasta "
	subprocess.run(combine_kmer_cmd, shell=True)

def blast_against_genomes():
	#create the combined genome database
	logging.debug("Creating genome database and BLASTing uniq k-mers.")
	combine_assemblies_cmd = f"cat assemblies/*fna > db.fna"
	subprocess.run(combine_assemblies_cmd, shell=True)
	create_db_cmd = f"makeblastdb -in db.fna -dbtype nucl"
	subprocess.run(create_db_cmd, shell=True)
	#blast uniq k-mers 
	blast_cmd = f"blastn -query tmp/diffKmers.fasta -db db.fna -outfmt 6 -qcov_hsp_perc 90 -perc_identity 90 -word_size 10 -out tmp/diff_kmer_blast_results.tsv -num_threads 30"
	subprocess.run(blast_cmd, shell=True)

def find_group_specific_kmers():
	logging.debug("Counting kmer presence in genome groups.")
	kmers = {}
	with open("tmp/diff_kmer_blast_results.tsv", 'r') as in_file:
		for line in in_file:
			line = line.rstrip().split("\t")
			if line[0] not in kmers:
				kmers[line[0]] = line[2]
				kmers[line[0]]['group1'] = 0
				kmers[line[0]]['group2'] = 0
			if group1_genomes[kmers[line[0]]] == 'g1':
				kmers[line[0]]['group1'] += 0
			if group2_genomes[kmers[line[0]]] == 'g2':
				kmers[line[0]]['group2'] += 0
			print(line[0])
			print(f"Group 1 count: kmers[line[0]]['group1']")
			print(f"Group 2 count: kmers[line[0]]['group2']")

def off_target_check(uniq_filtered_kmers, percent_identity):
	'''Use the kmer file filtered by frequency, GC content, melting temperature, and homopolymers as blast query against nt database.
	'''
	logging.debug(f"\tBLASTing selected k-mers against nt database to find off target hits.")
	blast_command = f"blastn -query tmp/{uniq_filtered_kmers} -db /storage/blastdb/v5/nt_v5 -outfmt '6 qaccver qseq saccver sscinames staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore' -qcov_hsp_perc {percent_identity} -perc_identity {percent_identity} -word_size 10 -out tmp/nt_blast.results -num_threads 20 -max_target_seqs 300"
	#print(blast_command)
	subprocess.run(blast_command, shell=True)

# nt_blast_results = "nt_blast.results"

def check_nt_blast_results(whiteList,blackList, nt_blast_results):
	'''
	Look at the k-mer to nt database blast results
	Use whitelist (bordetella and eventually accetable contaminants) and blacklist (bacteria and mammals) taxids
	Count hits to whitelist as on
	'''
	logging.debug(f"\tParsing off-target BLAST results selecting candidate primers.")
	kmer_blasthit_count = {}
	pass_taxids = [line.strip() for line in whiteList]
	no_pass_taxids = [line.strip() for line in blackList]
	with open(f"tmp/{nt_blast_results}", 'r') as in_file:
		for line in in_file:
			# print(line)
			line = line.rstrip().split('\t')
			if line[0] not in kmer_blasthit_count:
				kmer_blasthit_count[line[0]] = {}
				kmer_blasthit_count[line[0]]['pass_count'] = 0
				kmer_blasthit_count[line[0]]['no_pass_count'] = 0
			if line[4] in pass_taxids:
				kmer_blasthit_count[line[0]]['pass_count'] += 1
			elif line[4] in no_pass_taxids:
				kmer_blasthit_count[line[0]]['no_pass_count'] += 1
			else:
				kmer_blasthit_count[line[0]]['no_pass_count'] += 1
		percet_off_target = (int(kmer_blasthit_count[line[0]]['no_pass_count'])*100/int(kmer_blasthit_count[line[0]]['pass_count']))
		if percet_off_target < 10:
			with open("tmp/candidate_primers.fasta",'w') as out_file:
				out_file.write(f">{line[0]}\n{line[1]}")

	




# def check_dependencies():
# 	try:
# 		subprocess.run(["jellyfish"], stdout=devnull, stderr=devnull)
# 		print("Jellyfish: PASSED")
# 	except:
# 		sys.stderr.write("\n WARNING: Cannot find Jellyfish. Check that Jellyfish is downloaded and in your PATH\n\n")
# 		sys.exit()

# 	try:
# 		subprocess.run(['blastn'], stdout=devnull, stderr-devnull)
# 		print("Blast: PASSED")
# 	except:
# 		sys.stderr.write("\n WARNING: Cannot find Jellyfish. Check that Jellyfish is downloaded and in your PATH\n\n")
# 		sys.exit()




if __name__ == '__main__':
	main()


