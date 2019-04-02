#Borja :)
import sys
import numpy as np
from random import randint
args = sys.argv
input_file = None
output_file = None
def parse_args(args = None):
	global input_file, output_file
	if "-f" not in args:
		print "-f input_file (mandatory) -o output_file (optional) -m mutation_ratio -n num_copies -r read_size -c coverage"
		exit(1)
	input_file = args[args.index("-f")+1]
	output_file = args[args.index("-o")+1]

def change_reads():
	global input_file, output_file
	cont = 0
	if '1' in input_file.split('/')[-1]:
		extra = str('/1')
	else:
		extra = str('/2')
	with open(input_file, 'r') as f:
		with open(output_file, 'w+') as f2:
			for i,line in enumerate(f.readlines()):
				if i % 4 == 0:
					f2.write(line.split('\n')[0]+extra+'\n')
					if line.split('\n')[0][0]!='@':
						print line
				else:
					f2.write(line)

if __name__=='__main__':
	print 'Changing reads'
	parse_args(args)
	change_reads()
