import os
import sys
import time
import argparse
from collections import namedtuple
import shutil

try:
	from utils import utilitaire as ut, align, loop
except ModuleNotFoundError:
	import utilitaire as ut, align, loop


##
# Auxiliary function definitions
#

def processArgs():
	"""
	Returns argument parser for the main program.
	"""
	parser = argparse.ArgumentParser(description='Sequence Homology Iterative Finding Tool - version 1.5')
	## Mandatory arguments ###
	parser.add_argument('query', help='Query file in fasta format')
	parser.add_argument('target', help='Target file (search dataset) in fasta format')
	## General parameters ###
	parser.add_argument('-n', '--seqType', help='Sequence type (aa=default/nuc)', type=str, default='aa')
	parser.add_argument('-o', '--out_dir', help='Output directory (default="shift_result")', type=str, default='shift_result')
	parser.add_argument('-t', '--threads', help='Number of threads (default=1)', type=int, default=1)
	parser.add_argument('-r', '--runs', help='Maximum number of runs (default=1e6)', type=int, default=10**6)
	## Alignment parameters ###
	parser.add_argument('-a', '--aligner', help="Alignment program (blast=default/diamond/mmseqs)", type=str, default="blast")
	parser.add_argument('-p', '--pid', help='Identity threshold (in %%); default = 30.0, value > 0.0', type=float, default=30.0)
	parser.add_argument('-c', '--cover', help='Mutual coverage threshold (in %%); default = 80.0, value > 0.0', type=float, default=80.0)
	parser.add_argument('-e', '--e_value', help='E-value threshold; default = 1e-5, value > 0.0', type=float, default=0.00001)
	parser.add_argument('-m', '--mmseq_sensitivity', help='Sensitivity value (for mmseq2 only); default = 7.5, value > 0.0', type=float, default=7.5)
	## Iteration parameters ###
	parser.add_argument('-T', '--trim', help="Filter target by the min/max size of sequences in query (default True)", action="store_true", default=True)	# keep it or not ?
	## Module options ###
	parser.add_argument('-f', '--family', help='Check if input sequences form a family for the chosen parameters', action='store_true')
	parser.add_argument('-s', '--ssn', help='Compute sequences similarity network on all retrieved sequences (query + found)', action='store_true')
	## Output option ###
	parser.add_argument('-R', '--remove_aln', help="Remove aln directory and all its content; useful for limiting hard drive memory consumption", action='store_true')
	parser.add_argument('-K', '--keep_renum', help="Keep renumbered files and dictionaries", action='store_true')
	return parser


def pathSetup(args):
	"""
	Sets hardcoded (relative) paths for temporary and output files for the main program.
	Returns a namedtuple Path_tuple.
	"""
	isf_directory	   	= os.path.dirname(os.path.realpath(__file__))	   	# directory containing the current executable (here ihf.py) 		# what for ?
	working_dir	     	= os.path.abspath(args.out_dir)			 	# set by user: default cwd()/isf_wd
	db_dir		  	= os.path.join(working_dir, "db_dir")	       		# directory for storing databases (aligner)
	aln_out_dir	     	= os.path.join(working_dir, "aln_out_dir")	  	# directory for iterative alignment output files
	query_numID	     	= os.path.join(working_dir, "query_numID.faa")		# renumbered query file
	query_dico	     	= os.path.join(working_dir, "query_dico.txt")		# correspondence dictionary queryID -> query_numID
	filtered_target_numID   = os.path.join(working_dir, "filtered_target_numID.faa")# restricted and renumbered target file
	target_dico	     	= os.path.join(working_dir, "target_dico.txt")		# correspondence dictionary targetID -> target_numID
	target_db	       	= os.path.join(db_dir, "target_db")			# aligner target database
	log_file		= os.path.join(working_dir, "log_file.txt")		# log file containing most of the string output

	for d in [ working_dir, db_dir, aln_out_dir ]:
		ut.checkArg(d)

	Path_tuple = namedtuple("Path_tuple", ["working_dir", "db_dir", "aln_out_dir", "query_numID", "query_dico", 
					       "filtered_target_numID", "target_dico", "target_db", "log_file"])
		      
	myVars = Path_tuple(working_dir=working_dir, db_dir=db_dir, aln_out_dir=aln_out_dir, query_numID=query_numID, query_dico=query_dico, filtered_target_numID=filtered_target_numID,  
			    target_dico=target_dico, target_db=target_db, log_file=log_file)

	## Write arguments and parameters to the log_file ###
	to_print = ""
	this_format = "{: <50}\t{:<}\n"
	log = open(log_file, "a")
	to_print += this_format.format("executable", os.path.realpath(__file__))
	to_print += this_format.format("isf directory", os.path.abspath(isf_directory))
	to_print += this_format.format("working directory", os.path.abspath(working_dir))
	to_print += this_format.format("query", os.path.abspath(args.query))
	to_print += this_format.format("target", os.path.abspath(args.target))
	to_print += this_format.format("max number of runs", args.runs)
	to_print += "{: <50}\t{:<}%\n".format("pident threshold", args.pid)
	to_print += "{: <50}\t{:<}%\n".format("coverage threshold", args.cover)
	to_print += this_format.format("evalue threshold", args.e_value)
	to_print += "{: <50}\t{:<.1f}%\n".format("min size threshold",args.cover)
	to_print += "{: <50}\t{:<.1f}%\n".format("max size threshold",10**4/args.cover)
	to_print += this_format.format("aligner", args.aligner)
	to_print += this_format.format("seqType", args.seqType)
	to_print += "\n## Path variables\n"
	for name, value in myVars._asdict().items():
	       to_print += "    {: <50}\t{:<}\n".format(name, value)
	to_print += "\n"
	print(to_print)
	print(to_print, file=log)
	return myVars


def renameAlnFiles(seq_found, myVars):
	"""
	Recover original IDs for aln files -- using the query and target dico files -- and prepare fasta file of retrieved hits.
	"""
	print("Renaming alignment files...\n")
	ut.concat_file("large_dico.txt", myVars.query_dico, myVars.target_dico)						# aux file large_dico.txt contains both ID -> numID dictionaries
	aln_files_to_rename = [os.path.join(myVars.aln_out_dir, elem) for elem in os.listdir(myVars.aln_out_dir)]	# list of files to reID: alignment files
	for file in aln_files_to_rename:
		shutil.copy(file, file+'_num')										# keep aln files with numIDs
		ut.rename_blast(file, "large_dico.txt")									# transfer old ID back instead of the numID 
	print("Generating new query file...\n")
	ut.concat_file("large_numID.faa", myVars.query_numID, myVars.filtered_target_numID)				# join query and target in a unique aux large file
	new_query = os.path.join(myVars.working_dir, "shift_hits.faa")							# create new query file... (with shift_run:X appended after seq name)
	allSeq = dict()
	for key in ut.get_dico_num(myVars.query_dico):
		allSeq[key] = 0
	allSeq.update(seq_found)	   
	ut.grep_fasta_round(sequence_file = "large_numID.faa", round_found=allSeq, key_sort=lambda x:allSeq[x], out_file = new_query)
	allDict = ut.get_dico_num("large_dico.txt")
	seqDict = dict()
	for key in allSeq:
		seqDict[key] = allDict[key]+"\tshift_run:"+str(allSeq[key])
	ut.rename_fasta_from_dico(seqDict, new_query)									# ...finished new query file: transfer old IDs back
	ut.remove("large_numID.faa", "large_dico.txt")									# remove concatenated auxiliary files
	return new_query


##
# Main 
#

def main():
	parser = processArgs()
	args = parser.parse_args()
	## Test for alignment adequacy ###
	if args.aligner == 'diamond' and args.seqType == 'nuc':
	    print("Usage : DIAMOND only works with protein sequences -- try 'blast' or 'mmseqs' instead.")
	    sys.exit(1)
	## Set up files and directories ###
	query = os.path.abspath(args.query)			     	# absolute path of query file 
	target = os.path.abspath(args.target)			   	# absolute path of target file
	myVars = pathSetup(args)					# define paths of temporary hardcoded files
	os.chdir(myVars.working_dir)
	trim = args.trim
	seqType = args.seqType
	log = open(myVars.log_file, "a")
	## Start program ###
	starting_time = time.time()
	print("Starting SHIFT")
	### process_input sets all files and parameters for iterative step ======================================================================================
	db_size, min_size, max_size, q_index, t_index = ut.process_input(query=query, target=target, query_dico=myVars.query_dico, target_dico=myVars.target_dico, 
								target_out=myVars.filtered_target_numID, query_out=myVars.query_numID, 
								log_file=myVars.log_file, cov=args.cover, trim=trim)
	# process_input has computed : myVars.query_dico, myVars.target_dico, myVars.filtered_target_numID, myVars.query_numID ###
	### define aligner 
	aligner = align.Aligner(name=args.aligner, threads=args.threads, temp_dir=None, working_dir=myVars.working_dir, cover=args.cover, e_val=args.e_value, pid=args.pid, seqType=seqType, log=log)
	# =======================================================================================================================================================
	files_to_remove = []
	## Optional module 0 : gene family check ###
	if args.family:
		string_0 = "\nModule 0 : Checking query family composition -------------------\n"
		if q_index == 1:
			string_0 += "\nNothing to check: there is only one input sequence!\n"
		else:
			ssn_file = myVars.query_numID+".ssn"
			files_to_remove.append(ssn_file)
			### process_aln calls an all-vs-all alignment
			align.process_aln(query=myVars.query_numID, target=myVars.query_numID, aligner=aligner, out=ssn_file, out_dir=myVars.db_dir)
			# process_aln has computed : ssn_file
			CC = ut.computeCC(ssn_file)
			nb_CC = len(CC)
			if nb_CC > 1:
				 string_0 += "\n\tWarning: found %d connected components in the input SSN.\n\tThe input %d sequences in fasta file %s do not form a single gene family for these parameters.\n" % (nb_CC, q_index, query)
			else:
				 string_0 += "\n\tThe input %d sequences form a gene family for the specified parameters.\n" % q_index
		print(string_0)
		print(string_0, file=log)
		print("Module 0 finished -------------------\n")
	## Module 1 : iterative search ###
	string_1 = "Module 1 : iterative search -------------------\n"
	print(string_1)
	print(string_1, file=log)
	### make_db constructs the target database =========================================================================================================
	aligner.make_db(db_fasta=myVars.filtered_target_numID)	  
	# make_db computed : aligner-formatted target database in target_db
	### loop is the central function of the iterative module
	seq_found = loop.loop(query=myVars.query_numID, target=myVars.filtered_target_numID, aligner=aligner, aln_dir=myVars.aln_out_dir, run=args.runs)
	# loop computed : seq_found -- dict : numID -> round of iteration
	### remove all unnecessary files and db
	aligner.purge()
	#  =================================================================================================================================================
	print("Iteration done ++++++++++++ \n")
	### Rename and reorder data for output.
	if len(seq_found) > 0:
		if args.remove_aln:
			ut.remove(myVars.aln_out_dir)
		else:
			new_query = renameAlnFiles(seq_found, myVars)
	else:
		new_query = renameAlnFiles(seq_found, myVars)
	print("Module 1 finished ------------------- \n")
	## Optional module 2 : SSN of extended query ###
	if args.ssn:
		string_2 = "\nModule 2 : All-against-all alignment -------------------\n"
		print(string_2)
		all_against_all = os.path.join(myVars.working_dir, 'all_against_all.blastp')
		### process_aln calls all-vs-all alignment on new query
		align.process_aln(query=new_query, target=new_query, aligner=aligner, out=all_against_all, out_dir=myVars.db_dir, renum=True)
		#
		string_2_out = "\nAlignment output is {}".format(os.path.abspath(all_against_all))
		print(string_2_out)
		print(string_2+string_2_out, file=log)
		print("\nModule 2 finished ----------------- \n")
	## Ending ###
	total_time = time.time() - starting_time
	h, m, s = ut.time_parser(total_time)
	print("\t-------------------")
	print("SHIFT found: {} sequences".format(len(seq_found)))
	print("\t-------------------")
	outString = "Merged fasta of initial gene family and sequences found is {}".format(os.path.abspath(new_query))
	print(outString)
	print("\t-------------------")
	print("Fasta of initial seeds is {}".format(os.path.abspath(query)))
	# Cleaning up
	files_to_remove.extend([myVars.target_dico, myVars.filtered_target_numID, myVars.query_dico, myVars.query_numID, myVars.db_dir])
	if not args.keep_renum:
		ut.remove(*files_to_remove)
	print("\n"+outString,file=log)
	print("\nSHIFT took {} hours, {} minutes and {} seconds\n".format(str(h), str(m), str(s)), file=log)
	log.close()


if __name__ == "__main__":
	main()

