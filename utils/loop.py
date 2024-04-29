import os
import time
import argparse
from collections import defaultdict

try:
	from utils import utilitaire as ut, align, checkSegments
except ModuleNotFoundError:
	import utilitaire as ut, align, checkSegments


##
# Function definitions
#

def loop(query=None, target = None, aligner=None, aln_dir=None, run=None):
	"""
	Main program of the SHIFT software.
	Goal  : returns the sequences identified at successive iterations of alignment on the target dataset. 
	------------------
	Input : 
	query		: str		name of fasta input sequences (numerical IDs)
	target		: str		name of fasta target sequences
	aligner		: Aligner object
	aln_dir		: str		name of temporary directory for intermediate alignments
	run		: int		maximum number of iterations
	------------------
	Return : seq_found  dict(int:int)  dictionary : seq_ID -> run_number
	"""
	segments = defaultdict(set)
	seq_found = dict()
	run_number = 1
	cont = 1
	current_query = query												# start value of current query (fasta file)
	while run_number <= run and cont:
		current_aln = os.path.join(aln_dir, '_' + str(run_number))						# path of output file: aln_out_dir/_X
		aligner.align(fasta = current_query, output = current_aln, ban_list = set(seq_found.keys()))
		segments = checkSegments.select_new_sequences(clean_aln = current_aln, threshold = aligner.cover/100, run_number = run_number, segments = segments)
		current_query = 'ihf_temp.fa'
		ut.grep_fasta_num(sequence_file = target, id_ = segments.keys(), out_file = current_query) 		# current_query: update file
		for seq in segments:											# found sequences: update dictionary 
			seq_found[seq] = run_number
		runString = 'run_number: {}, new sequences found: {}, total sequences found: {}'.format(run_number, len(segments), len(seq_found))
		print(runString)
		print(runString, file=aligner.log)
		run_number += 1
		cont = len(segments)
	if os.stat(current_aln).st_size == 0:
		ut.remove(current_aln)
	if run_number >= 2:				# avoid blunder if current_query is still equal to query...
		ut.remove(current_query)
	return seq_found

##
# Standalone usage procedures
#

def processArgs():
	"""
	Returns parser for loop standalone usage.
	"""
	parser = argparse.ArgumentParser(description='SHIFT - version 1.0')
	parser.add_argument('query', help='Fasta file of query gene families')
	parser.add_argument('target', help='Fasta file of target homolog sequences')
	parser.add_argument('-aligner', help="Define alignment program (blast=default/diamond/mmseq2)", type=str, default="blast")
	parser.add_argument('-th', help='Number of threads to use', type=int, default=1)
	parser.add_argument('-out_dir',	help="Output directory")	
	parser.add_argument('-pid', help='Threshold limit for identity; default = 30.0, value > 0.0', type=float, default=30.0)
	parser.add_argument('-cov', help='Threshold limit for mutual coverage; default = 80.0, value > 0.0', type=float, default=80.0)
	parser.add_argument('-e_val', help='Threshold limit for e-value; default = 1e-5, value > 0.0', type=float, default=0.00001)

	return parser

def main():
	"""main function for standalone use -- not fully functional yet?"""
	parser = processArgs()
	args = parser.parse_args()
	query = args.query
	if args.out_dir:
		this_dir = os.path.join(os.getcwd(), args.out_dir)
		os.makedirs(this_dir, exist_ok = True)
	else:
		this_dir = os.getcwd()
	alignerName = args.aligner
	aligner = align.Aligner(name=args.aligner, threads=args.th, cover=args.cov, e_val=args.e_val, pid=args.pid)
	aligner.make_db(args.target)
	loop(query=query, target=args.target, aligner=aligner, aln_dir=this_dir, run=1e6)

if __name__ == "__main__":
	main()

