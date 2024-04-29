import os
import re
from copy import deepcopy
from collections import defaultdict
import operator
import shutil


##
# Functions: sequence segments 
#


def check_overlap(aln, qA, qB):				
	"""
	Central routine of segments compatibility condition : compute intersection of alignment match from aln with segment (qA, qB) on query sequence
	Input : aln 	str (blastp outfmt 6 format line)
	      : qA, qB  int
	Return: sA, sB  int (projection of intersection segment (qA, qB) with (sstart, send) onto search sequence)
		overlap int (length of the intersection)
	"""
	qstart, qend, sstart, send = int(aln[5]), int(aln[6]), int(aln[7]), int(aln[8])
	qseq, sseq = aln[12:]
	aln_len = len(qseq)
	match = {}
	q, s = qstart, sstart
	dA = 0
	flag = True
	for i in range(aln_len):
		lq = qseq[i]
		ls = sseq[i]
		dq, ds = (lq != '-'), (ls != '-')
		if (dq + ds) == 2:	# if there is a match
			match[q] = s
			if q <= qB:		# shift the end while it is smaller than qB
				end = q
			if q >= qA and flag:	# start is the first match index that is larger than qA
				start = q
				flag = False
		q += dq
		s += ds
		if flag:
			dA += 1
	q1, q2, s1, s2 = start, end, match[start], match[end]
	overlap = q2 - q1 + 1
	return s1, s2, overlap


def is_included(segm1, segm2):
	"""
	Input : segm = (l_min, start, end)
	Returns 0 if and only if segm2 is strictly included in segm1, 1 if they are equal and 2 otherwise.
	H: len(segm1) >= len(segm2)
	"""
	L1, a1, b1 = segm1
	L2, a2, b2 = segm2
	if a1 == a2 and b1 == b2:
		return 1
	elif a1 <= a2 and b2 <= b1:
		return 0
	else:
		return 2


def keep_maximal(segment_set):
	"""
	Input: set_of_segments
	segment : S = (l_min, start, end)
	maximality condition : S' > S iff (S included in S') or (S = S' and S'.l_min <= S.l_min)
	uses : is_included for trichotomy condition
	Returns : subset of maximal segments from segment_set
	"""
	return_list = []
	# Sort list by length (decreasing) then by l_min value (ascending)
	sorted_segments = sorted(list(segment_set), key = lambda x:(x[2]-x[1]+1, -x[0]))
	while sorted_segments:
		segm2 = sorted_segments.pop()
		incl = 2
		ct = 0
		M = len(return_list)
		while incl >= 2 and ct < M:
			segm1 = return_list[ct]
			incl = is_included(segm1, segm2)
			ct += 1
		if incl == 2:
			return_list.append(segm2)
	return set(return_list)


def refine(segments):
	"""
	Input : dictionary seq -> segment_set
	Returns : dictionary seq -> subset_of_maximal_segments of segment_set
	"""
	newSegments = dict()
	for seq in segments:
		newSegments[seq] = keep_maximal(segments[seq])
	return newSegments


##
# Functions: loop program
#

def select_new_sequences(clean_aln, threshold, run_number, segments):
	"""
	Wrapper around the 3 options for manage_sequences.
	--------------------
	Goal  : selects the target sequences from clean_aln that satisfy the compatibility condition with respect to the current segments dictionary, 
	Returns : updated segments dictionary. Target sequences are contained in newSegments.keys()
	-------------------
	Input :
	clean_aln : a blastp-formatted (outfmt 6) alignment file  qseqid sseqid evalue pident bitscore qstart qend sstart send qlen slen
						          	     0     1       2     3       4       5      6    7      8    9   10
	threshold : float
	run_number : int
	segments : dict(set_of_tuples)
	Returns :
	newSegments : dict(set_of_tuples).
	"""
	new_aln = clean_aln+".temp"
	newSegments = manage_sequences(clean_aln, new_aln, threshold, run_number, segments)
	shutil.move(new_aln, clean_aln)
	return newSegments


def manage_sequences(clean_aln, new_aln, threshold, run_number, segments=None):
	"""
	mode "all". 
	segments dictionary structure : seq_ID -> set((len_min, start, end)) : set of segments with value len_min and end points (start, end) over seq_ID.
	----------	
	Transitive threshold condition : 
	a match (A, B, startA, endA, startB, endB) from clean_aln is kept
	if there exists a segment (len_min, start, end) in segments[A] such that the intersection of (start, end) and (startA, endA) 
	1) covers at least threshold % (default 80%) of sequence A.
	2) is at least len_min
	-> item (l_min, sB, eB) in newSegments[B] where:
	* l_min = max( len_min, threshold * length(B) )
	* (sB, eB) = projection of intersection on B -- return value of check_overlap 
	"""
	newSegments = defaultdict(set)
	if run_number == 1:							# select all targets from clean_aln at the first round 
		with open(clean_aln) as inp, open(new_aln, "w") as out:
			for line in inp:
				aln = line.strip().split()
				query, search, sstart, send, qlen, slen = int(aln[0]), int(aln[1]), int(aln[7]), int(aln[8]), int(aln[9]), int(aln[10])
				l_min = threshold * max(qlen, slen)
				newSegments[search].add((l_min, sstart, send))	# take all segments (including the overlapping ones)
				out.write(line)
	else:									# at next rounds, filter by compatibility condition
		with open(clean_aln) as inp, open(new_aln, "w") as out:
			for line in inp:
				aln = line.strip().split()
				query, search, sstart, send, qlen, slen = int(aln[0]), int(aln[1]), int(aln[7]), int(aln[8]), int(aln[9]), int(aln[10])

				segment_list = list(segments[query])		# list of all matching segments (l_min, qstart, qend) of query

				cpt = 0
				flag = False					
				while cpt < len(segment_list):					# loop over segments on query
					len_min, qA, qB = segment_list[cpt]			# len_min = transitive length	qA = match start in query	qB = match end in query 
					sA, sB, overlap = check_overlap(aln, qA, qB)		# sA, sB = candidates for projection bounds on the search seq
					cpt += 1
					if overlap/qlen >= threshold and overlap >= len_min:	# threshold condition and transitive condition
						flag = True
						l_min = max(len_min, threshold * slen)		# update transitive length value
						newSegments[search].add((l_min, sA, sB))	# include all segments passing threshold and transitive conditions (including the overlapping ones)
				if flag:
					out.write(line)						# write aln only once
	newSegments = refine(newSegments)							# keep maximal segments only
	return newSegments
