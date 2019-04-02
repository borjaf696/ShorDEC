#Borja :)
import numpy as np
import sys

def simulate_new_haplotype(fasta,new_fasta, divergence):
	ids = []
	sequences = []
	with open(fasta,'r') as f:
		for line in f.readlines():
			if line[0] == '>':
				ids.append(line.strip())
			else:
				sequences.append(line)
	new_sequences = []
	nt = ['A','C','G','T']
	for sequence in sequences:
		new_sequence = ''
		for p in sequence:
			new_sequence+=(p if np.random.binomial(1,divergence,1) == 0 else nt[np.random.randint(0,3,1)[0]])
		new_sequences.append(new_sequence)
	with open(new_fasta,'w+') as f:
		for i,n_sequence in enumerate(new_sequences):
			f.write(ids[i]+'_'+str(divergence)+'\n')
			f.write(n_sequence+'\n') if i < len(new_sequences)-1 else f.write(n_sequence)
		
if __name__=="__main__":
	print 'Creating new mutant haplotypes'
	fasta = sys.argv[1]
	divergence = float(sys.argv[2])
	new_fasta = ('/').join(fasta.split('/')[:-1])+'/'+fasta.split('/')[-1].split('.')[0]+'_'+str(divergence)+'.'+fasta.split('.')[-1]
	print 'Fasta: ', fasta
	print 'NewFile: ',new_fasta 
	simulate_new_haplotype(fasta, new_fasta, divergence)
