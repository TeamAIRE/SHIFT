import subprocess
import argparse
import os
import re
import sys
import math
import time
import shutil
import tempfile
from collections import defaultdict

try:
	from utils import utilitaire as ut
except ModuleNotFoundError:
	import utilitaire as ut


##
# Auxiliary function definitions 
#

def present(ident, banList):
	'''T * list[T] -> bool
	Returns true if ident is present in banList.
	Designed to circumvent the issue of identifier's type in banList'''
	try:
		res = int(ident) in banList
	except ValueError:
		res = ident in banList
	return res


def filterAln(aln, cover, e_val, pid, ban_list):
	'''str * float**3 * list[T] -> bool
	Returns True if aln has to be filtered out.
	Version for diamond or mmseqs aligners.''' 
	return (not aln or present(aln[1], ban_list) or aln[0] == aln[1])


def filterBlast(aln, cover, e_val, pid, ban_list):
	'''str * float**3 * list[T] -> bool
	Returns True if aln has to be filtered out.
	Version for blastp.''' 
	tcov = 100*(int(aln[8]) - int(aln[7]) + 1)/int(aln[10])
	qcov = 100*(int(aln[6]) - int(aln[5]) + 1)/int(aln[9])
	return (not aln or present(aln[1], ban_list) or aln[0] == aln[1] or tcov < cover or qcov < cover or float(aln[3]) < pid or float(aln[2]) > e_val)


def myWhich(name):
	'''str -> str
	Returns path to binary corresponding to the command "name" or exits if non-existing.
	'''
	binary = shutil.which(name)
	if not binary:
		print("Cannot find executable %s: install it, or check your PATH" % name)	
		sys.exit(1)
	return binary


##
# Class definitions

class Aligner:
	"""
	This class implements an object allowing a convenient handling of alignment software. 
	It can currently be used with blastp, diamond and mmseqs.
	"""

	def __init__(self, name=None, threads=None, temp_dir=None, working_dir=None, cover=0, e_val=10, pid=0, mm=None, db_fasta=None, seqType="aa", log=None):
		"""
		Constructor for class Aligner.
		name		:	blast/diamond/mmseq2
		threads		:	nb of threads for parallel aligner run
		temp_dir	:	optional temporary directory (for mmseq2)
		working_dir	:	directory (aligner database storage)
		cover		:	mutual coverage threshold
		e_val		:	maximum e_value threshold
		pid		:	minimum % id threshold
		mm		:	sensitivity parameter (mmseqs only)
		db_fasta	:	optional fasta file to construct aligner database.
		seqType 	:	sequence type (aa/nuc) -> sets binary for blast and option for mmseqs; incompatible with diamond
		"""
		# Set path to aligner binary
		if name == 'blast':
			if seqType == "aa":
				name = 'blastp'
			else:
				name = 'blastn'
		if name == 'diamond' and seqType == 'nuc':
			print("DIAMOND only works with protein sequences -- try 'blast' or 'mmseqs' instead")
			sys.exit(1)
		binary = myWhich(name)
		# Set working directory
		if not working_dir:
			working_dir = os.getcwd()			
		working_dir = os.path.join(working_dir,'db_dir')
		self.log = log
		# additional parameters required by mmseqs only
		if name == 'mmseqs':
			if mm:
				mm = float(mm)
			else:
				print('Using mmseqs sensitivity value 7.5; set other value (with option -m) if desired, and run again\n')
				mm = 7.5
			temp_dir = os.path.join(working_dir, "mmseq_temp_dir")
		elif name != 'blastp' and name !='blastn' and name != 'diamond':
			print(" -- Wrong type or name of aligner")
			sys.exit(1)
		# Set the attributes
		self.name = name
		self.binary = binary
		self.threads = threads

		self.temp_dir = temp_dir 
		self.sensitivity = mm
		self.working_dir = working_dir

		self.pid = pid
		self.e_val = e_val
		self.cover = cover
		self.mx = 100000
		self.seqType = seqType

		self.ban_list = set([])
	
		if db_fasta:
			self.db = self.make_db(db_fasta)
		

	def make_db(self, db_fasta, alt=False):
		"""
		Formats a fasta file as a database out_db in the format suitable for self.
		"""
		t_name = os.path.basename(db_fasta+"_db")
		out_db = os.path.join(self.working_dir, t_name)
		if self.name == 'blastp':
			mkdb = myWhich("makeblastdb")
			cmd = '{} -dbtype prot -hash_index -parse_seqids -in {} -out {}'.format(mkdb, db_fasta, out_db)
		elif self.name == 'blastn':
			mkdb = myWhich("makeblastdb")
			cmd = '{} -dbtype nucl -hash_index -parse_seqids -in {} -out {}'.format(mkdb, db_fasta, out_db)
		elif self.name == 'diamond':
			cmd = "{} makedb --in {} --db {}".format(self.binary, db_fasta, out_db)
			out_db += ".dmnd"
		elif self.name == 'mmseqs':
			cmd = '{} createdb {} {}'.format(self.binary, db_fasta, out_db, self.threads)
		else:
			print(" ** Wrong type of aligner")
			sys.exit(1)		
		dbString1 = "Making {} database with the command line:\n\t {}".format(self.name, cmd)
		dbString2 = "Created database with radical : {}\n".format(out_db)
		child = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
		child.wait()
		print(dbString1)
		print(dbString2)
		if not alt:
			# default behaviour: sets the aligner.db attribute
			self.db = out_db
		else:
			# alternative : for mmseqs, since it requires also a db for the query ...
			return out_db


	def align(self, fasta, output, db=None, max_seq_target=None, cover=None, e_val=None, pid=None, ban_list=None, directed=True):
		'''Performs the alignment of fasta with aligner and produces file "output".
		The alignment is free from self-hits and filtered by cover, e_value, pid and best reciprocal hit (when directed is False).
		By default, uses the attributes of the aligner, which can be individually overriden.'''
		if not db:
			db = self.db
		if not max_seq_target:
			max_seq_target = self.mx
		if not cover:
			cover = self.cover
		if not e_val:
			e_val = self.e_val
		if not pid:
			pid = self.pid
		if not ban_list:
			ban_list = self.ban_list
		else:
			self.ban_list = ban_list
		number_of_sequences = sum(1 for line in open(fasta) if line.startswith(">"))
		nbString = 'There are {} sequences to align with {}\n'.format(number_of_sequences, self.name)

		if self.name == 'blastp':
			cmd = ['{}'.format(self.binary), '-query', '{}'.format(fasta), '-db', '{}'.format(os.path.abspath(db)),
			   '-out', '{}'.format(output), '-evalue', '{}'.format(e_val), '-num_threads', '{}'.format(self.threads),
			   '-outfmt', "6 qseqid sseqid evalue pident bitscore qstart qend sstart send qlen slen gaps qseq sseq",
		        "-seg", "yes", "-soft_masking", "true", "-max_target_seqs", "10000"]
			child = subprocess.Popen(cmd, shell=False, stdout=subprocess.DEVNULL)
			child.wait()
			
		elif self.name == 'blastn':
			cmd = ['{}'.format(self.binary), '-query', '{}'.format(fasta), '-db', '{}'.format(os.path.abspath(db)),
			   '-out', '{}'.format(output), '-evalue', '{}'.format(e_val), '-num_threads', '{}'.format(self.threads),
			   '-outfmt', "6 qseqid sseqid evalue pident bitscore qstart qend sstart send qlen slen gaps qseq sseq",
		        "-soft_masking", "true", "-max_target_seqs", "10000"]
			child = subprocess.Popen(cmd, shell=False)#, stdout=subprocess.DEVNULL)
			child.wait()

		elif self.name == 'diamond':
			cmd = "{} blastp --db {} --query {} --out {} --id {} --outfmt 6 qseqid sseqid evalue pident bitscore qstart qend sstart send qlen slen gaps qseq_gapped sseq_gapped --threads {} \
--query-cover {} --subject-cover {} --evalue {} --max-target-seqs {}".format(self.binary, db, fasta, output, pid, self.threads, cover, cover, e_val, max_seq_target)
			child = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			child.wait()

		elif self.name == 'mmseqs':
			# alignment with mmseq needs the query db as a database
			with tempfile.TemporaryDirectory(prefix='temp_split_dir', dir=self.working_dir) as tmp_db_dir:
				db_fasta = os.path.join(tmp_db_dir, 'mseq_db_query')
				aln_db = os.path.join(tmp_db_dir, '.aln_db')
				db_fasta = self.make_db(fasta, alt=True)
				if self.seqType == "aa":
					cmd1 = "{} search {} {} {} {} -s {} --min-seq-id {} -c {} -e {} --max-seqs 100000 --max-rejected 100 --threads {} --alignment-mode 3 -a".format(self.binary,
							db_fasta, db,  aln_db, self.temp_dir, self.sensitivity, pid/100 , cover/100, e_val, self.threads)
				else:
					cmd1 = "{} search {} {} {} {} -s {} --min-seq-id {} -c {} -e {} --max-seqs 100000 --max-rejected 100 --threads {} --alignment-mode 3 --search-type  3 -a".format(self.binary, db_fasta, db,  aln_db, self.temp_dir, self.sensitivity, pid/100 , cover/100, e_val, self.threads)

				child = subprocess.Popen(cmd1, shell=True, stdout=subprocess.DEVNULL)
				child.wait()
				cmd2 = "{} convertalis {} {} {} {}  --threads {} --format-output \"query,target,evalue,pident,bits,qstart,qend,tstart,tend,qlen,tlen,gapopen,qaln,taln\"".format(
					self.binary, db_fasta, db, aln_db, output, self.threads)
				child = subprocess.Popen(cmd2, shell=True, stdout=subprocess.DEVNULL)
				cmd = cmd1+"\n"+cmd2
				child.wait()
		else:
			print(" $$ Wrong type of aligner")
			sys.exit(1)
		alString = "Launching {} with the following command:\n{}\n".format(self.name, cmd)
		print(alString)		# only on the stdout
		print(nbString)
		self.clean(output, directed = directed)


	def clean(self, align_file, directed = True):
		'''Submethod of align. Removes self-hits and filters by mutual cover, e_value and pid, as well as best reciprocal e_value hit if directed is False.'''
		temp_file = align_file + ".temp_aln_file"
		id_0 = None
		dico_best_hit = {}
		if self.name == 'blastp' or self.name == 'blastn':
			filt = filterBlast
		elif self.name == 'diamond' or self.name == 'mmseqs':
			filt = filterAln
		else:
			print("Wrong type of aligner")
			sys.exit(1)
		print("Cleaning",align_file)
		with open(align_file) as in_f, open(temp_file, "w") as out:
			for line in in_f:
				aln = line.strip().split("\t")
				if not aln or filt(aln, self.cover, self.e_val, self.pid, self.ban_list):
					continue
				query, search = aln[0:2]
				if id_0 != query:
					id_0 = query
					for k, v in dico_best_hit.items():
						out.write("{}\n".format("\t".join(v)))
					dico_best_hit = {(query, search): aln}
				else:
					key = (query, search)
					if key not in dico_best_hit:
						dico_best_hit[key] = aln
					else:
						old_aln = dico_best_hit[key]
						if float(old_aln[2]) > float(aln[2]):  # compare e-values
							dico_best_hit[key] = aln
			if dico_best_hit:
				for k, v in dico_best_hit.items():
					out.write("{}\n".format("\t".join(v)))
		if not directed:
			temp_file = process_undirected_aln(temp_file)
		shutil.move(temp_file, align_file)
	
	def purge(self):
		'''Remove all db files of aligner.'''
		print("Removing database files",self.db)
		if self.name == 'diamond':
			ut.remove(self.db)
		elif self.name == 'blastp' or self.name == 'mmseqs':
			for filename in os.listdir(self.working_dir):
				if re.search(self.db, os.path.abspath(filename)):
					ut.remove(filename)



##
# Principal functions
#

def process_undirected_aln(align_file):
	"""
	Select best e_value hit for all reciprocal hits in align_file (intended for all-against-all align files). 
	"""
	temp_file = align_file + ".temp"

	id_0 = None
	dico_best_hit = {}
	already_seen = defaultdict(int)
	num_lines = sum(1 for line in open(align_file))
	ct = int(round(math.sqrt(num_lines)))				# dump lines after every batch of ct 

	with open(align_file) as in_f, open(temp_file,"w") as dump:
		count = 0
		for line in in_f:
			count += 1
			spt = line.strip().split("\t")
			if not spt or spt[0] == spt[1]:			# ignores empty lines or self-loops : in principle, is useless now...
				continue
			head, tail = min(spt[0], spt[1]), max(spt[0], spt[1])			
			if id_0 != spt[0]:						# reference sequence: all following lines match this sequence until it changes
				if id_0:
					already_seen[id_0] = 1					# id_0 as query cannot appear any longer
				if count % ct == 0:
					toDel = set()
					for k, v in dico_best_hit.items():			# when ref changes, we dump the content of dico_best_hit into the outFile...
						if already_seen[k[0]] and already_seen[k[1]]:
							dump.write("{}\n".format("\t".join(v)))
							toDel.add(k)
					for k in toDel:
						del dico_best_hit[k]
				id_0 = spt[0]
				dico_best_hit[(head, tail)] = [str(head), str(tail)] + spt[2:]
			else:
				key = (head, tail)				# if ref is the same, we create a key (ref,match) ...
				if key not in dico_best_hit:				# ... and either store the result if it's the only one so far...
					dico_best_hit[key] = [str(head), str(tail)] + spt[2:]
				else:							# ... or compare it according to the e-value
					old_spt = dico_best_hit[key]
					if float(old_spt[2]) < float(spt[2]):			## New version: passes "best" filter
						dico_best_hit[key] = [str(head), str(tail)]+spt[2:]
		if dico_best_hit:							# ... and dump the last matching in the outFile
			for k, v in dico_best_hit.items():
				dump.write("{}\n".format("\t".join(v)))
	shutil.move(temp_file, align_file)
	return align_file
	

def process_aln(query = None, target = None, aligner = None, out = None, out_dir = None, cover = None, e_val = None, pid = None, renum = False):
	'''Main procedure of ALIGN.
	query	:	fasta input
	target	:	fasta optional input (otherwise = query)
	aligner	:	str (blastp/diamond/mmseqs)
	out	:	str (output file name)
	out_dir	:	str (output dir)
	cover	:	float (minimum reciprocal cover)
	e_val	:	float (maximum evalue)
	pid	:	float (minimum % identity)
	renum	:	bool (switches identifiers with integers and back)
	'''
	if renum:
		print("Option renum on")
		ut.make_dictionary_and_get_size(in_file = query, out_file = query+"_num", out_dico = query+"_dico")
		ut.make_dictionary_and_get_size(in_file = target, out_file = target+"_num", out_dico = target+"_dico")
		the_query, the_target = query+"_num", target+"_num"
	else:
		print("Option renum off")
		the_query, the_target = query, target
	t_name = os.path.basename(the_target+"_db")
	the_db = os.path.join(out_dir, t_name)
	print("query = %s\ntarget = %s\ndb = %s" % (the_query, the_target, the_db))
	aligner.make_db(db_fasta = the_target)
	#print(aligner.db)
	outfile = os.path.join(out_dir, out)
	aligner.align(fasta = the_query, output = outfile, directed = (query != target))
	to_remove = [the_db]
	if renum:
		if query != target:
			ut.concat_file(query+"_dico", target+"_dico")
		ut.rename_blast(in_file = outfile, dico_file = query+"_dico")
		to_remove.extend([query+"_num", query+"_dico", target+"_num", target+"_dico"])
	for the_file in to_remove:
		print("Removing file %s" %  the_file)
		ut.remove(the_file)
	aligner.purge()
	return 0


##
# Standalone use procedures
#


def processArgs():
	parser = argparse.ArgumentParser(description='Sequence Homolog Iterative Finding Tool version 1.0')
	parser.add_argument('query', help='Fasta file of query gene families')
	parser.add_argument('-target', help='Fasta file of target gene families')
	parser.add_argument('-aligner', help="Define alignment program (blastp=default/diamond/mmseqs)", type=str, default="blastp")
	parser.add_argument('-out', help="Outfile name")
	parser.add_argument('-out_dir',	help="Output directory")	
	parser.add_argument('-th', help='Number of threads to use', type=int, default=1)
	parser.add_argument('-pid', help='Threshold limit for identity; default = 30.0, value > 0.0', type=float, default=30.0)
	parser.add_argument('-cov', help='Threshold limit for coverage; default = 80.0, value > 0.0', type=float, default=80.0)
	parser.add_argument('-e_val', help='Threshold limit for e-value; default = 0.00001, value > 0.0', type=float, default=0.00001)
	parser.add_argument('-renum', help="Renumber query and target fasta files", action='store_true')
	parser.add_argument('-nuc', help="Nucleotide", action='store_true')
	return parser


def main():
	'''Reads the arguments, sets the aligner object and launches the procedures.'''
	parser = processArgs()
	args = parser.parse_args()
	query = args.query
	if args.target:
		target = args.target
	else:
		target = args.query
	if args.out:
		out = args.out
	else:
		out = "_".join([query, args.aligner, "aln"])
	if args.out_dir:
		this_dir = os.path.join(os.getcwd(), args.out_dir)
		os.makedirs(this_dir, exist_ok = True) 
	else:
		this_dir = os.getcwd()
	if args.nuc:
		seqType = "nuc"
	else:
		seqType = "aa"
	print("Initial arguments: ", args, "\nquery=%s, target=%s, out=%s, out_dir=%s" % (query, target, out, this_dir))
	aligner = Aligner(name=args.aligner, threads=args.th, temp_dir=None, working_dir=this_dir, cover=args.cov, e_val=args.e_val, pid=args.pid, seqType=seqType)
	start_time = time.time()
	process_aln(query=query, target=target, aligner=aligner, out=out, out_dir=this_dir,  renum=args.renum)
	print("Alignment output is {}".format(out))
	total_time = time.time() - start_time
	print("Alignment performed in {:.2f} s.".format(total_time))
	
if __name__ == '__main__':
	main()

