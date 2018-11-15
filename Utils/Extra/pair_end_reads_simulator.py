#Borja :)
import sys
import numpy as np
from random import randint
args = sys.argv

mutation_ratio = 0.1
fail_on_reads = 0.05
len_var = 0.01
num_copies = 4
output_file = "output.fasta"
input_file = None
read_size = 150
insert_size = 100
max_delta = 10
coverage = 25

nuc = ['A', 'C', 'G', 'T']

def parse_args(args = None):
	global input_file, mutation_ratio, output_file, num_copies, read_size, coverage
	if "-f" not in args:
		print "-f input_file (mandatory) -o output_file (optional, with out .fasta) -m mutation_ratio -n num_haplotypes -r read_size -c coverage -i insert_size -d max_delta -g read_len_var"
		exit(1)
	if "-m" in args:
		mutation_ratio = float(args[args.index("-m")+1])
	if "-n" in args:
		num_copies = int(args[args.index("-n")+1])
	if "-o" in args:
		output_file = args[args.index("-o")+1]
		if (".fasta" in output_file):
			output_file = output_file[0:len(output_file-7)]
	if "-r" in args:
		read_size = int(args[args.index("-r")+1])
	if "-c" in args:
		coverage = int(args[args.index("-c")+1])
	if "-i" in args:
		insert_size = int(args[args.index("-i")+1])
	if "-d" in args:
		max_delta = int(args[args.index("-d")+1])
	if "-g" in args:
		len_var = float(args[args.index("-g")+1])
	input_file = args[args.index("-f")+1]

def _new_read_copies(read):	
	global mutation_ratio, num_copies
	reads = []
	for i in range(num_copies):
		l_read = list(read[0:len(read)-1])
		l_read = [l_read[i] if t==0 else nuc[randint(0,3)] for i,t in enumerate(np.random.binomial(1,mutation_ratio,len(l_read)))]
		reads.append(("").join(l_read))
	return reads

def new_reads(reads):
	#Creating haplotypes
	haplotypes = []
	for read in reads:
		haplotypes = haplotypes + _new_read_copies(read)
	return haplotypes

def _add_fails_read(read):
	global fail_on_reads
	l_read = list(read[0:len(read)-1])
	l_read = [l_read[i] if t==0 else nuc[randint(0,3)] for i,t in enumerate(np.random.binomial(1,fail_on_reads,len(l_read)))]
	read = ("").join(l_read)

def _cut_read(read):
	global read_size, insert_size, max_delta, len_var
	left_reads, right_reads = [], []
	left_read, right_read = "",""
	for i in range(0,len(read)-insert_size-read_size):
		insert_local = insert_size+randint(-max_delta, max_delta)
		diff = read_size+randint(int(-read_size*len_var), int(read_size*len_var))
		last_len = min(i+diff, len(read))
		i_left = min(last_len+insert_local, len(read))
		left_bound = min(i_left+diff, len(read))
		left_read = read[i:last_len]
		right_read = read[i_left:left_bound]
		if (len(left_read) > len(right_read)):
			for i in range(len(right_read),len(left_read)):
				right_read  = right_read + 'N'
		left_reads.append(left_read)
		right_reads.append(right_read)
		left_read, right_read = "",""
	return left_reads, right_reads

def add_errors(reads):
	global coverage
	final_reads = []
	for read in reads:
		final_reads = final_reads + [read]*coverage
	print(len(final_reads))
	for read in final_reads:
		_add_fails_read(read)
	right_reads, left_reads = [], []
	for read in final_reads:
		left, right = _cut_read(read)
		left_reads = left_reads + left
		right_reads = right_reads + right
	print(len(left_reads)," ",len(right_reads))
	return left_reads, right_reads

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
	reads_left, reads_right = add_errors(reads)
	with open(output_file+"_1.fasta",'w+') as f:
		for i,read in enumerate(reads_left):
			f.write(">"+str(i)+"\n")
			f.write(read+"\n")
	with open(output_file+"_2.fasta","w+") as f:
		for i, read in enumerate(reads_right):
			f.write(">"+str(i)+"\n")
			f.write(read+"\n")
	print "END!"
