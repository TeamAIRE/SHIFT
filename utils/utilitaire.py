import os
import sys
import subprocess
import shutil
import math
import igraph


def concat_file(dest, *args):
	"""
	Concatenate the list of files in args as dest.
	Appends the result if dest already exists.
	"""
	with open(dest, 'ab') as wfd:
		for f in args:
			with open(f, 'rb') as fd:
				shutil.copyfileobj(fd, wfd, 1024*1024*10)
				wfd.write(b"\n")


def rename_blast(in_file, dico_file):
	"""
	Overwrites blast output in_file with original sequence names as recorded in dico_file.
	H : Keys and values in the dico file are treated as strings
	"""
	dico = get_dico_from_dico_file(dico_file)
	temp_file = in_file + ".tmp"
	with open(in_file) as in_f, open(temp_file, 'w') as out:
		for line in in_f:
			spt = line.strip().split()
			spt[0] = dico[spt[0]]
			spt[1] = dico[spt[1]]
			out.write("{}\n".format("\t".join(spt)))
	shutil.move(temp_file, in_file)


def get_dico_from_dico_file(dico_file):
	"""
	:param dico_file:	"\t" formatted several column file where identifier is the last field in each row : string"\t"identifier
	:return:		dictionary identifier -> string 
	"""
	dico = {}
	with open(dico_file) as inp:
		for line in inp:
			spt = line.strip().split('\t')
			dico[spt[-1]] = '  '.join(spt[0:-1])
	return dico


def get_dico_num(dico_file):
	"""
	Alternative version of previous function with integer keys.
	"""
	dico = {}
	with open(dico_file) as inp:
		for line in inp:
			spt = line.strip()
			if spt:
				spt = spt.split('\t')
				dico[int(spt[-1])] = '  '.join(spt[0:-1])
	return dico


def grep_fasta(sequence_file, id_, out_file):
	"""
	if sequence have similar name in the fasta it may not get them all
	this function will retrieve n sequence then return n being the number of uniq id
	:param sequence_file: a fasta file
	:param id_: can be a file one id line or an iterable
	:param out_file: output a fasta file
	:return: None
	"""

	set_id = set()
	if isinstance(id_, str):
		if os.path.isfile(id_):
			with open(id_) as id_1:
				for line in id_1:
					set_id.add(line.strip())
	else:
		set_id = id_

	nb_seq = len(set_id)
	cpt = 0
	with open(sequence_file, 'r') as seqfile, open(out_file, 'w') as out_:
		for prompt, sequ in get_seq_one_by_one(seqfile):
			if prompt in set_id:

				out_.write(">{}\n{}\n".format(prompt, chunck_sequence(sequ)))
				cpt += 1
			if cpt == nb_seq:
				return

def grep_fasta_num(sequence_file, id_, out_file):
	"""
	if sequence have similar name in the fasta it may not get them all
	this function will retrieve n sequence then return n being the number of uniq id
	:param sequence_file: a fasta file
	:param id_: can be a file one id line or an iterable : that is rather weird (!)
	:param out_file: outputs a fasta file
	:return: None
	"""

	set_id = set()
	if isinstance(id_, str):
		if os.path.isfile(id_):
			with open(id_) as id_1:
				for line in id_1:
					set_id.add(line.strip())
	else:
		set_id = id_

	nb_seq = len(set_id)
	cpt = 0
	with open(sequence_file, 'r') as seqfile, open(out_file, 'w') as out_:
		for prompt, sequ in get_seq_one_by_one(seqfile):
			if int(prompt) in set_id:

				out_.write(">{}\n{}\n".format(prompt, chunck_sequence(sequ)))
				cpt += 1
			if cpt == nb_seq:
				return

def grep_fasta_round(sequence_file, round_found, key_sort, out_file):
	"""
	if sequence have similar name in the fasta it may not get them all
	this function will retrieve n sequence then return n being the number of uniq id
	:param sequence_file: a fasta file
	:param id_: can be a file one id line or an iterable
	:param out_file: outputs a fasta file
	:return: None
	"""

	outDict = dict()

	with open(sequence_file, 'r') as seqfile:
		for prompt, sequ in get_seq_one_by_one(seqfile):
			if int(prompt) in round_found:
				outDict[int(prompt)] = sequ

	with open(out_file, 'w') as out_:
		for key in [k for k in sorted(round_found, key=key_sort)]:
			prompt = str(key)
			sequ = outDict[key]
			out_.write(">{}\n{}\n".format(prompt, chunck_sequence(sequ)))


def get_seq_one_by_one(in_file):
	"""
	Generator returning lists [name, sequence] for each sequence in in_file
	"""
	sequence = ''
	prompt = ''
	for line in in_file:
		if line.startswith('>'):
			if sequence:
				yield [prompt, sequence]
			sequence = ''
			prompt = line.strip()[1:]
		else:
			sequence += line.strip()
	yield [prompt, sequence]


def chunck_sequence(sequence):
	inc = 80
	cpt = 0
	seq = ""
	while cpt <= len(sequence):
		seq += sequence[cpt:cpt + inc] + '\n'
		cpt += inc
	return seq.strip()


def rename_fasta(dico_file, fasta):
	"""
	Rewrites fasta file with original names contained in dico_file
	Warning : overwrites original fasta file
	"""
	my_dico = get_dico_from_dico_file(dico_file)
	temp_fasta = fasta + ".temp"
	with open(temp_fasta, 'w') as my_temp_file, open(fasta) as fasta_:
		for prompt, seq in get_seq_one_by_one(fasta_):
			my_temp_file.write(">{}\n{}\n".format(my_dico[prompt], chunck_sequence(seq)))
	shutil.move(temp_fasta, fasta)


def rename_fasta_from_dico(dico, fasta):
	"""
	Rewrites fasta file with original names contained in dico_file
	Warning : overwrites original fasta file
	"""
	temp_fasta = fasta + ".temp"
	with open(temp_fasta, 'w') as my_temp_file, open(fasta) as fasta_:
		for prompt, seq in get_seq_one_by_one(fasta_):
			my_temp_file.write(">{}\n{}\n".format(dico[int(prompt)], seq))
	shutil.move(temp_fasta, fasta)

def remove(*args):
	for element in args:
		try:
			if os.path.isfile(element):
				os.remove(os.path.abspath(element))
			elif os.path.isdir(element):
				shutil.rmtree(os.path.abspath(element))
		except:
			pass


def split_fasta(fasta, tmp_dir, n_split, fasta_split):
	fasta_split_ok = True
	try:
		child = subprocess.Popen(["{}".format(fasta_split), '-f',  '{}'.format(fasta), "-o", "{}".format(tmp_dir), "-c",  "{}".format(n_split)])
		child.wait()
	except:
		print("w")
		fasta_split_ok = False

	if not fasta_split_ok:
		print("fastasplit seems not present splitting file will take longer")
		liste_handler = []
		for indice in range(n_split):
			liste_handler.append(open(os.path.join(tmp_dir, "{}.faa".format(indice)), 'w'))
		generator = loop_counter(n_split)
		with open(fasta) as f_in:
			for prompt, seq in get_seq_one_by_one(f_in):
				indice = generator.__next__()
				liste_handler[indice].write(">{}\n{}\n".format(prompt, chunck_sequence(seq)))
		for elem in liste_handler:
			elem.close()


def loop_counter(end, start=0):
	count = start
	while True:
		yield count
		count += 1
		if count >= end:
			count = start


def process_input(query, target, query_dico, target_dico, target_out, query_out, log_file, cov, trim=None):
	"""
	Processes query and target files. Returns number of target sequences retained and their min/max lengths, as well as the boundary indices of renumbering of sequences.
	Outputs query_dico, target_dico, target_out and query_out files.
	"""
	# first make dictionary of the query and return the min max len of query sequences
	min_len, max_len, indice_dico = make_dictionary_and_get_size(in_file = query, out_file = query_out, out_dico = query_dico)
	min_sequence_len = int(math.ceil(min_len * (cov / 100))) 
	max_sequence_len = int(math.floor(max_len * (100 / cov)))
	# filter target by query length
	cpt, target_index, last_index = screen_target(min_sequence_len = min_sequence_len, max_sequence_len = max_sequence_len,	indice_dico = indice_dico, 
							in_file = target, out_dico = target_dico, out_fasta = target_out, trim = trim)
	to_print = ""
	to_print += "Gene Family/Query:\n"
	to_print += "\tShortest sequence: {} aa\n\tLongest sequence: {} aa\n".format(min_len, max_len)
	to_print += "\tFirst sequence index: {}\n\tLast sequence index: {}\n".format(1, indice_dico-1)
	to_print += "Out of the {} sequences from the target input:\n".format(target_index)
	to_print += "{} lie in the length range of the query sequences and are hereafter selected as the target database\n".format(cpt)
	if trim:
		to_print += "\tMin length kept: {}\n\tMax length kept: {}\n".format(min_sequence_len, max_sequence_len)
	to_print += "\tFirst sequence index: {}\n\tLast sequence index: {}\n".format(indice_dico, last_index-1)
	print(to_print)
	add_to_file(log_file, to_print)
	if cpt <= 0 :
		print("No sequences selected")
		raise AssertionError
	return (cpt, min_len, max_len, indice_dico-1, last_index-1)


def make_dictionary_and_get_size(in_file, out_file, out_dico, start = 1):
	"""
	Renumber fasta in_file as out_file and store correspondence in out_dico file. Returns min and max lengths of sequences and last index value.
	"""
	first_bool = True
	id_indice = start  # id start for one by default
	with open(in_file, buffering=1000000) as input_file, open(out_file, 'w') as output_file, open(out_dico, 'w') as dico_file:
		for prompt, sequence in get_seq_one_by_one(input_file):
			prompt = prompt.replace("\t", "  ")
			# This part aims to find the min max sequence length
			length_seq = len(sequence)
			if first_bool:
				min_seq = length_seq
				max_seq = length_seq
				first_bool = False
			if length_seq < min_seq:
				min_seq = length_seq
			if length_seq > max_seq:
				max_seq = length_seq
			# Create the dictionary
			dico_file.write("{}\t{}\n".format(prompt, id_indice))
			output_file.write(">{}\n{}\n".format(id_indice, chunck_sequence(sequence)))
			# Increment the value of the indice
			id_indice += 1
	return (min_seq, max_seq, id_indice)


def screen_target(min_sequence_len, max_sequence_len, indice_dico, in_file, out_dico, out_fasta, trim):
	"""
	Filter fasta in_file by min_ and max_sequence_len, and renumber as out_file with correspondence dictionary out_dico. 
	Returns an interger triplet: nb_of_retained_sequences, nb_of_input_sequences, last_index_of_renumbering.
	"""
	cpt = 0
	nb = 0
	with open(in_file, buffering=10000000) as input_file, open(out_dico, 'w') as dico_file, open(out_fasta, 'w') as fasta_out:
		for prompt, sequence in get_seq_one_by_one(input_file):
			# if the sequence is in the threshold we write
			if trim:
				if min_sequence_len <= len(sequence) <= max_sequence_len:
					dico_file.write("{}\t{}\n".format(prompt, indice_dico))
					fasta_out.write(">{}\n{}\n".format(indice_dico, chunck_sequence(sequence)))
					cpt += 1
					indice_dico += 1
			else:
				dico_file.write("{}\t{}\n".format(prompt, indice_dico))
				fasta_out.write(">{}\n{}\n".format(indice_dico, chunck_sequence(sequence)))
				cpt += 1
				indice_dico += 1
			nb += 1
	return cpt, nb, indice_dico


def add_to_file(in_file, txt):
	with open(in_file, "a") as out:
		out.write("{}\n\n".format(txt))



def output_round(found, query_dic, target_dic, round_file):
	"""
	Write found dictionary : seq_numID -> run_number (from main_loop.isf) as round_file with original identifiers (as stored in query_dic and target_dic) 
	"""
	q_dico = get_dico_from_dico_file(query_dic)
	t_dico = get_dico_num(target_dic)					# this function reads the dictionary "backwards"
	with open(round_file, "w") as f:
		f.write("id\tround\n")							# header
		for item in q_dico:							# query sequence + '0'
			f.write("{}\t{}\n".format(q_dico[item], '0'))			
		for elem in found:							# target seq + 'run'
			f.write("{}\t{}\n".format(t_dico[elem], str(found[elem])))	


def time_parser(second):
	f_hour, f_minute, f_sec = 0,0,0
	f_minute = int(second / 60.0)
	f_sec = second % 60
	f_hour = int(f_minute / 60)
	f_minute = f_minute % 60

	return (f_hour, f_minute, round(f_sec,2))


def checkArg(arg):
	try:
		os.mkdir(arg)
	except FileExistsError:
		print("File/directory %s exists: check your arguments or remove dir/file" % arg)
		sys.exit(1)


def computeCC(edgeFile):
	tmp = edgeFile+".tmp"
	cmd = """cut -f1-2 %s > %s""" % (edgeFile,tmp)
	os.system(cmd)
	g = igraph.read(tmp, format="ncol", directed=False, names=True)
	comp = g.clusters(mode='weak')
	remove(tmp)
	return comp


## Debug

def examineFile(filename, k=10):
	nb = 0
	f = open(filename)
	while nb <= k:
		line = f.readline().strip()
		print(line)
		nb += 1
	f.close()

def examineFileAll(filename):
	f = open(filename)
	for line in f:
		line = f.readline().strip()
		print(line)
	f.close()

