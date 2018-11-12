#Borja :)
import sys
import numpy as np
from random import randint
args = sys.argv

mutation_ratio = 0.1
fail_on_reads = 0.01
num_copies = 5
output_file = "output.fasta"
input_file = None
read_size = 150
coverage = 100

nuc = ['A', 'C', 'G', 'T']

def parse_args(args = None):
	global input_file, mutation_ratio, output_file, num_copies, read_size, coverage
	if "-f" not in args:
		print "-f input_file (mandatory) -o output_file (optional) -m mutation_ratio -n num_copies -r read_size -c coverage"
		exit(1)
	if "-m" in args:
		mutation_ratio = int(args[args.index("-m")+1])
	if "-n" in args:
		num_copies = int(args[args.index("-n")+1])
	if "-o" in args:
		output_file = args[args.index("-o")+1]
	if "-r" in args:
		read_size = int(args[args.index("-r")+1])
	if "-c" in args:
		coverage = int(args[args.index("-c")+1])
	input_file = args[args.index("-f")+1]

def _fail_reads(read):
	#Min_coverage = coverage
	global read_size, coverage
	reads = []
	for i in range(coverage):
		for j in range(len(read)-read_size+1):
			l_read = list(read[j:j+read_size])
			l_read = [l_read[i] if t == 0 else nuc[randint(0,3)] for i,t in enumerate(np.random.binomial(1,fail_on_reads,len(l_read)))]
			reads.append(("").join(l_read))
	return reads	
	

def add_errors(reads):
	final_reads = []
	for read in reads:
		final_reads = final_reads + _fail_reads(read)
	return final_reads

def _new_read_copies(read):	
	global mutation_ratio, num_copies
	reads = []
	for i in range(num_copies):
		l_read = list(read[0:len(read)-1])
		l_read = [l_read[i] if t==0 else nuc[randint(0,3)] for i,t in enumerate(np.random.binomial(1,mutation_ratio,len(l_read)))]
		reads.append(("").join(l_read))
	return reads

def new_reads(reads):
	final_reads = []
	for read in reads:
		final_reads = final_reads + _new_read_copies(read)
	return final_reads

if __name__ == "__main__":
	parse_args(args)
	reads = []
	with open(input_file,'r') as f:
		for line in f.readlines():
			if ">" not in line:
				reads.append(line)
	reads = new_reads(reads)
	parse = output_file.split('/')
	parse[-1] = "haplotypes_"+parse[-1]
	parse = "/".join(parse)
        with open(parse,'w+') as f:
		for i, read in enumerate(reads):
			f.write("Haplotype number: "+str(i)+"\n")
			f.write(read+"\n")
	reads =	add_errors(reads)
	with open(output_file,'w+') as f:
		for i,read in enumerate(reads):
			f.write(">"+str(i)+"\n")
			f.write(read+"\n")
	print "END!"
	
	
	






