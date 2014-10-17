"""
pull out all subsequence aligning against a user-specified reference
"""

import sys, os, re
from time import time
import hyphyAlign

# tried out the Biopython pairwise alignment module - it's a memory hog :-P
#from Bio import pairwise2

gap_prefix = re.compile('^[-]+')
gap_suffix = re.compile('[-]+$')


def get_boundaries (str):
	# return a tuple giving indices of subsequence without gap prefix and suffix
	res = [0,len(str)]
	left = gap_prefix.findall(str)
	right = gap_suffix.findall(str)
	if left:
		res[0] = len(left[0])
	if right:
		res[1] = len(str) - len(right[0])
		
	return res


def convert_fasta (lines):	
	blocks = []
	sequence = ''
	for i in lines:
		if i[0] == '$': # skip h info
			continue
		elif i[0] == '>' or i[0] == '#':
			if len(sequence) > 0:
				blocks.append([h,sequence])
				sequence = ''	# reset containers
				h = i.strip('\n')[1:]
			else:
				h = i.strip('\n')[1:]
		else:
			sequence += i.strip('\n').upper()
	blocks.append([h,sequence])	# handle last entry
	return blocks


def main():
	if len(sys.argv) < 3:
		print 'Usage: python parse_fasta.py [FASTA] [ref seq] [outfile]'
		sys.exit()
	
	no_genome = False
	if '-nogenome' in sys.argv:
		no_genome = True
	
	# import sequences from FASTA file
	fn = sys.argv[1]
	infile = open(fn, 'rU')
	fasta = convert_fasta(infile.readlines())
	infile.close()
	
	print 'read', len(fasta), 'sequences from', fn
	
	"""
	# find the longest sequence
	max_slen = 0
	which_max = None
	for i in range(len(fasta)):
		h, s = fasta[i]
		slen = sum( [s.count(nuc) for nuc in 'ACGT'] )
		if slen > max_slen:
			max_slen = slen
			which_max = i
	
	long_h, long_s = fasta[which_max]
	print 'longest sequence %s has %d chars' % (long_h, max_slen)
	"""
	
	# alignment settings
	
	hyphy = hyphyAlign.HyPhy._THyPhy (os.getcwd(), 1) # instance of HyPhy
	dump = hyphy.ExecuteBF('MESSAGE_LOGGING = 0;', False)
	
	hyphyAlign.change_settings(hyphy, 
								alphabet = hyphyAlign.nucAlphabet, 
								scoreMatrix = hyphyAlign.nucScoreMatrix, 
								gapOpen = 20,
								gapOpen2 = 20,
								gapExtend = 10,
								gapExtend2 = 10,
								noTerminalPenalty = 1)
	
	
	refseq = sys.argv[2].upper()
	outfile = open(sys.argv[3], 'w')
	
	for h, s in fasta:
		if no_genome and 'genome' in h:
			continue
		
		if 'Modified' in h or 'Patent' in h:# or len(s) > 10000:
			print 'skipping ' + h
			continue
		
		query, ref, score = hyphyAlign.pair_align(hyphy, refseq, s)
		left, right = get_boundaries(ref)
		if score > len(refseq):
			print h
			outfile.write('>'+h+'\n'+query[left:right]+'\n')
	
	outfile.close()
	
	
	
	"""
	# do quick-and-dirty pairwise alignment against longest sequence
	start_time = time()
	
	outfile = open(sys.argv[2], 'w')
	outfile.write('>'+long_h+'\n'+long_s+'\n')
	
	for h, s in fasta:
		if h == long_h:
			continue
		print h
		query, ref, score = hyphyAlign.pair_align(hyphy, long_s, s)
		#alignments = pairwise2.align.localxx(long_s, s)
		#query, ref = alignments[0][:2]
		
		match = gap_prefix.findall(query)
		if match:
			left = len(match[0])
		
		outfile.write('>'+h+'\n' + '-'*left + s + '\n')
	
	outfile.close()
	
	end_time = time()
	print 'that took %d seconds' % (end_time - start_time)
	"""
	

if __name__ == "__main__":
	main()
