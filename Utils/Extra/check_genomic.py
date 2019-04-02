#/usr/bin 
import sys
def load_truth(path):
	set_truth = []
	with open(path, 'r') as f:
		for line in f.readlines():
			set_truth.append(line.strip())
	return set_truth

def get_freqs(path, truth_set):
	list_freqs =[0]*20000
	with open(path,'r') as f:
		for line in f.readlines():
			line = line.strip()
			line_bucket = line.split()
			if int(line_bucket[1]) == 1 or int(line_bucket[1]) == 2:
				continue
			if line_bucket[0] in truth_set:
				list_freqs[int(line_bucket[1])] += 1
	return list_freqs

def to_csv(freqs):
	with open('genomic.csv','w+') as f:
		for i,freq in enumerate(freqs):	
			f.write(str(i)+';'+str(freq)+'\n')

if __name__=="__main__":
	print 'lets do this'
	path_truth = sys.argv[1]
	path_check = sys.argv[2]
	print 'Path: ', path_truth,' ',path_check
	truth_set = load_truth(path_truth)
	list_freqs = get_freqs(path_check, truth_set)
	to_csv(list_freqs)
