#!/bin/usr

truth = '../../solids_truth'
trial = '../../solids'

if __name__=="__main__":
	print 'Lets work'
	realkmers, falsekmers = [],[]
	with open(truth, 'r') as f:
		for l in f:
			realkmers.append(l)
	print 'RealKmers: ', len(realkmers)
	hit,total = 0,0
	with open(trial, 'r') as f:
		for l in f:
			total += 1
			if l in realkmers:
				hit += 1
				realkmers.remove(l)
			else:
				falsekmers.append(l)
	print 'Hits: ', hit,' Total: ',total
	print 'Sensibility: ', hit/len(realkmers)
	print 'Precission: ',hit/total
