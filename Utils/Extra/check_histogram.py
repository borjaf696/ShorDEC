#/usr/bin 
import sys

def get_freqs(path):
	count = 0
	freqs = []
	with open(path, 'r') as f:
		for line in f.readlines():
			line = line.strip()
			line = line[:-1] if line[-1] == ',' else line
			if str.isdigit(line):
				line = int(line)
				if count % 2 ==1:
					freqs.append(line)
				count += 1
	return freqs

def to_csv(freqs):
	with open('histogram.csv','w+') as f:
		for i,freq in enumerate(freqs):	
			f.write(str(i)+';'+str(freq)+'\n')

if __name__=="__main__":
	print 'lets do this'
	path = sys.argv[1]
	print 'Path: ', path
	freqs = get_freqs(path)
	to_csv(freqs)
