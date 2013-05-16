#Copyright (C) 2013 Robert Lanfear
#
#This program is free software: you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by the
#Free Software Foundation, either version 3 of the License, or (at your
#option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#General Public License for more details. You should have received a copy
#of the GNU General Public License along with this program.  If not, see
#<http://www.gnu.org/licenses/>.


def combinations(iterable, r):
	pool = tuple(iterable)
	n = len(pool)
	if r > n:
		return
	indices = range(r)
	yield tuple(pool[i] for i in indices)
	while True:
		for i in reversed(range(r)):
			if indices[i] != i + n - r:
				break
		else:
			return
		indices[i] += 1
		for j in range(i+1, r):
			indices[j] = indices[j-1] + 1
		yield tuple(pool[i] for i in indices)


def import_fasta_as_dict(alignment_filename):
	######################        Import fasta file as dict        #################################
	alignment_file = open(alignment_filename, 'r')
	#and set up the alignment as a dictionary
	seqdict = {}
	#get the seqs from the alignment file
	start = 0 # a cheap trick to do this simply
	for line in alignment_file:
		if line.startswith(">"): #it's a new sequence
			#first, stash the old sequence into the dictionary
			if start == 1:
				seqdict[seq_name] = seq
			seq_name = line
			seq = '' #an empty sequence
		else: #we must have a sequence, maybe with a newline on the end of it
			seq_part = line.strip("\n") #remove newlines as we go
			seq = ''.join([seq, seq_part])
		start=1
	#add the last sequence, which you'll have missed because I'm too lazy to write this properly
	seqdict[seq_name] = seq	
	alignment_file.close()
	
	return seqdict

def get_protein_domains(sequence):
	'''
	###################################################################################################
	################################## Calculate TM Domains ###########################################
	#
	#Take a string of an amino acid sequence, return a list of Amino Acid objects
	#Each object has 2 things: the amino acid, and what it's predicted to be:
	#	TM: transmembrane domain
	#	EX: extracellular domain
	#	IN: intracellular domain
	#	NA: anything else (like signal peptides)
	#
	#This function uses the phobius web server, here: http://phobius.sbc.su.se/index.html 
	#
	#If you use it, you should cite the appropriate refrences, listed here: http://phobius.sbc.su.se/instructions.html
	#
	'''

	class AminoAcid():
		def __init__(self, seq, region):
			self.seq = seq	#which amino acid it is
			self.region = region #is it a TM, EX, or IN

	from ClientForm import ParseResponse
	from urllib2 import urlopen
		
	#go visit the website and get the result
	response = urlopen("http://phobius.binf.ku.dk/index.html")
	forms = ParseResponse(response, backwards_compat=False)
	form = forms[0]
	form["protseq"] = sequence
	result = urlopen(form.click())
	
	#process that result to a string of amino acids
	amino_list = []
	counter = 0
	for line in result.readlines():
		#first figure out what kind of structure we're talking about
		words = line.split()
		if len(words)>0 and words[0] == "FT":
			start = int(words[2]) - 1 #zero indexed
			end = int(words[3]) - 1 #zero indexed
			I_am = "NA"
			if words[1] == "TRANSMEM": 
				I_am = "TM"
			if words[1] == "TOPO_DOM":
				if words[4] == "CYTOPLASMIC.":
					I_am = "IN"
				if words[4] == "NON":
					I_am = "EX"
			
			if words[1] == "TRANSMEM" or words[1] == "SIGNAL" or words[1] == "TOPO_DOM":
				r = end - start + 1
				for i in range(r):
					aa = AminoAcid(sequence[counter], I_am)
					amino_list.append(aa)
					counter += 1
	
	return amino_list

def write_reduced_mase_alignment(alignment, list_of_tips_to_keep):
	
	''' A dumb function which takes as input an alignment object (from AlignIO)
	and a list of tip names, and outputs a mase alignment called "aln.mase", including only 
	those species which were included in the list of tip names. 
	
	I really haven't tried very hard here, but put this in just in case I ever needed it again
	'''
	tips = list_of_tips_to_keep
	
	#make a dict from the nexus alignment
	seqs = {}
	for thing in alignment._records:
		seqs[thing.id] = thing
	
	# make sub-alignment
	globalheader=";; created by Rob's python function reduce_alignment_mase()"
	localheader=";no comment"
	sub_alignment = [globalheader]
	for tip in tips:
		seq = seqs[tip].seq.data
		line = '\n'.join([localheader, tip, seq])
		sub_alignment.append(line)

	#write file
	alignoutfile = open('aln.mase', 'w')
	
	num_tax = len(tips) #this might at least guard against some more obvious potential errors, where the number of species in the list of tips for whatever reason doesn't match what ends up in the alignment... let's hope.
	length = len(seq)
	
	for line in sub_alignment:
		alignoutfile.write('%s\n' %(line))
	alignoutfile.close()
	
	print "Done, %d taxa with sequence length %d written to file 'aln.mase'" %(num_tax, length)

	return seqs

def reduce_alignment(alignment_name, list_of_tips_to_keep, format):
	
	''' A dumb function which takes as input a filename for an alignment, a format for the alignment
	and a list of tip names, and outputs a phylip alignment of the same name, including only 
	those species which were included in the list of tip names. 
	
	I really haven't tried very hard here, but put this in just in case I ever needed it again
	'''
	
	from Bio import AlignIO
	
	tips = list_of_tips_to_keep
	
	alignment = AlignIO.read(open(alignment_name, "rU"), format)
	
	#make a dict from the nexus alignment
	seqs = {}
	for thing in alignment._records:
		seqs[thing.id] = thing
	
	# make sub-alignment
	sub_alignment = []
	for tip in tips:
		seq = seqs[tip].seq.data
		line = '   '.join([tip, seq])
		sub_alignment.append(line)

	#write file
	filename = alignment_name.split('.')[0]
	filename = '.'.join([filename, 'phy'])
	alignoutfile = open(filename, 'w')
	
	num_tax = len(tips) #this might at least guard against some more obvious potential errors, where the number of species in the list of tips for whatever reason doesn't match what ends up in the alignment... let's hope.
	length = len(seq)
	
	alignoutfile.write('%d %d\n' %(num_tax, length))
	
	for line in sub_alignment:
		alignoutfile.write('%s\n' %(line))
	alignoutfile.close()
	
	print "Done, %d taxa with sequence length %d written to file: %s" %(num_tax, length, filename)

	return seqs
	
def multiply_brlens_by_random_variable(tree, distribution, parameters):

	'''
	Take a given newick tree, with branch lengths of some sort
	this is designed to multiply trees by rates from a given distribution
	
	The function will multiply each branchlength in the tree by a branchlength from
	a given distribution.
	
	Input: phylip tree as a string; and a distribution and parameters as follows:
		'betavariate' [alpha, beta]
		'expovariate' [lambd]
		'gammavariate' [alpha, beta]
		'gauss' [mu, sigma]
		'lognormvariate' [mu, sigma]
		'normalvariate' [mu, sigma]
		'paretovariate' [alpha]
		'uniform' [min, max]
		
	'''
	import random

	def _get_random_number(dist, params):	
		if dist == 'betavariate':
			num = random.betavariate(params[0], params[1])
		if dist == 'expovariate':
			num = random.expovariate(params[0])
		if dist == 'gammavariate':
			num = random.gammavariate(params[0], params[1])
		if dist == 'gauss':
			num = random.gauss(params[0])
		if dist == 'lognormvariate':
			num = random.lognormvariate(params[0], params[1])
		if dist == 'normalvariate':
			num = random.normalvariate(params[0], params[1])
		if dist == 'paretovariate':
			num = random.paretovariate(params[0])
		if dist == 'uniform':
			num = random.uniform(params[0], params[1])
		return num

	tree_as_list = []
	for thing in tree:
		tree_as_list.append(thing)
	i = 0
	
	new_tree = []
	
	while i<len(tree_as_list):
		if tree_as_list[i] == ':':
			new_tree.append(tree_as_list[i])
			i+=1	
			#and you know that what follows a colon is a number...
			num = []
			while tree_as_list[i].isdigit() or tree_as_list[i] =='.':
				num.append(tree_as_list[i])
				i+=1
			num = float(''.join(num))
			num = str(float(num*_get_random_number(distribution, parameters)))
			new_tree.append(num)
		else:
			new_tree.append(tree_as_list[i])
			i+=1
	new_tree = ''.join(new_tree)
	return new_tree
	
	

def assign_random_brlens(tree, distribution, parameters):
	
	'''
	Take a given phylip tree as input, and assign random branchlengths

	The tree can't be labelled or anything like that
	
	Input a phylip tree and the name of the distribution and any variables it takes as a list:
		'betavariate' [alpha, beta]
		'expovariate' [lambd]
		'gammavariate' [alpha, beta]
		'gauss' [mu, sigma]
		'lognormvariate' [mu, sigma]
		'normalvariate' [mu, sigma]
		'paretovariate' [alpha]
		'uniform' [min, max]
		
	'''
	import random

	def _get_random_number(dist, params):	
		if dist == 'betavariate':
			num = random.betavariate(params[0], params[1])
		if dist == 'expovariate':
			num = random.expovariate(params[0])
		if dist == 'gammavariate':
			num = random.gammavariate(params[0], params[1])
		if dist == 'gauss':
			num = random.gauss(params[0])
		if dist == 'lognormvariate':
			num = random.lognormvariate(params[0], params[1])
		if dist == 'normalvariate':
			num = random.normalvariate(params[0], params[1])
		if dist == 'paretovariate':
			num = random.paretovariate(params[0])
		if dist == 'uniform':
			num = random.uniform(params[0], params[1])
		return num
		
	tree_as_list = []
	
	for thing in tree:
		if thing == ',':
			thing = ":%f," %(_get_random_number(distribution, parameters))
		if thing == ')':
			thing = ":%f)" %(_get_random_number(distribution, parameters))
	
		tree_as_list.append(thing)
		
	newtree = ''.join(tree_as_list)
	
	return newtree
	
def generate_labelling_permutation(tree, localclocks):
	
	'''takes a phylip tree and a list of local clocks and returns a tree with branches randomly labelled
		
	the list of local clocks should be the number of branches in each, e.g. [10, 20] would be
	two local clocks one with 10 branches in, and the other with 20 in
	
	N.B. The remaining branches in the tree constitute a third local clock, but theres no need
	to specify this (and if you try to, it will mess things up)
	
	
	'''
	from random import shuffle
	
	#first make a new tree as a list, and a list of branches
	#in the new tree, the each branch will be labelled like this '##branch_x'
	#where x is an integer. Hopefully this won't conflict with any taxon labels...
	newtree = []
	branches = []
	
	branchnumber = 1
	
	for thing in tree:
		
		branch = ''
				
		if thing == ',' or thing == ')':
			#it's a branch
			branch = "%s%d" %('##branch_', branchnumber)
			branches.append(branch)
			branchnumber += 1
		
		newtree.append(branch)
		newtree.append(thing)
	
	#print newtree
	#print branches
	
	#now randomise the list of branches, twice just to be sure
	shuffle(branches)
	shuffle(branches)

	#make a dict of branches with labels they'll have
	branch_num = 0
	branchdict = {}
	for i in range(len(localclocks)):		
		label = '#%d' %(i+1)
		#label the appropriate number of branches for each number in localclocks list
		for j in range(localclocks[i]):
			branchdict[branches[branch_num]] = label
			branch_num += 1

	#now replace all branches in the tree with the appropriate label
	#if the branch isn't in branchdict, that means it doesn't get a lable (and PAML will assign it to rate category zero)
	for thing in branches:
		if branchdict.has_key(thing):
			newtree[newtree.index(thing)] = branchdict[thing]
		else:
			newtree.remove(thing)
	
	#now turn the newtree list into a phylip looking tree
	outputtree = ''.join(newtree)
	
	#print outputtree
	
	return outputtree

def permute_branch_labels(tree):	
	'''
	takes a phylip tree with branch labels, and shuffles them around	
	all unlabelled branches will be left unlabelled
	N.B. only works with single digit labels, e.g. #0 to #9,
	all labelling categories are shuffled amongst each other.

	if using for the LCPT, the 'nuisance rate' category should remain unlabelled
	that way it doesn't get shuffled, but PAML will still assign those branches their
	own rate category.

	'''
	from random import shuffle

	
	#first make a new tree as a list, and a list of branches
	#in the new tree, the each branch will be labelled like this '##branch_x'
	#where x is an integer. Hopefully this won't conflict with any taxon labels...
	newtree = []
	branches = []
	# a dict of branch labels, keyed by branch_id
	branch_labels = {}
	is_labelled = 0
	branchnumber = 1
	i = 0
	for i in range(len(tree)):
		thing = tree[i]
		#reset the branch, just in case
		branch_id = ''
		if thing == '#':
			#flag the branch labels as they appear, this flag
			#gets reset to zero after you find the next branch
			is_labelled = 1
			#extract the label - this will ONLY work for single-digit labels
			branch_label = '#%s' %(tree[i+1])
		if thing == ',' or thing == ')':
			#it's a branch
			#the list of branch labels is added to the new tree, just to demarcate which branches are which
			branch_id = "%s%d" %('##branch_', branchnumber)
			branches.append(branch_id)
			if is_labelled == 1:
				branch_labels[branchnumber] = branch_label
			branchnumber += 1
			is_labelled = 0
		newtree.append(branch_id)
		if is_labelled == 0:
			#in this way, the new tree is free of all labels
			newtree.append(thing)
	#now get a list of the labelled branches, and shuffle it around
	label_list = branch_labels.keys()
	shuffle(label_list)
	shuffle(label_list)
	#make new dict of labelled branches with the shuffled list - equivalent to shuffling labels among branches
	new_branch_labels = {}
	i = 0
	for key in branch_labels:
		new_branch_labels[branches[key-1]] = branch_labels[label_list[i]]
		i += 1
	#now replace all branches in the tree with the appropriate label
	#if the branch isn't in branchdict, that means it doesn't get a lable (and PAML will assign it to rate category zero)
	for thing in branches:
		if new_branch_labels.has_key(thing):
			newtree[newtree.index(thing)] = new_branch_labels[thing]
		else:
			newtree.remove(thing)
	#now turn the newtree list into a phylip looking tree
	outputtree = ''.join(newtree)
	return outputtree

def bootstrap(list):
	'''
	bootstrap a list of numbers to create a new list of the same length
	'''
	
	import random
	
	newlist = []
	
	for i in range(len(list)):
		newlist.append(random.choice(list))
		
	return newlist

def bootstrap_coeff_of_var(list, reps):
	"""
	calculate the coefficient of variation, and the 95% CIs by bootstrapping 'reps' times
	"""
	
	import sys
	sys.path.append( "/Users/Rob/Documents/python_scripts/" )
	from RobPhyloFunctions import bootstrap
	from numpy import mean, std
	
	obs_cov = std(list)/mean(list)

	if reps<100:
		print "Less than 100 replications is a bit pointless, so resetting to do 100 replications"
		reps = 100

	upperlimit = int(round(reps*0.975))
	lowerlimit = int(round(reps*0.025))
	
	covs = []
	
	for i in range(reps):
		newlist = bootstrap(list)
		covs.append(std(newlist)/mean(newlist))
		
	covs.sort()
	
	upperCI = covs[upperlimit]
	lowerCI = covs[lowerlimit]
	
	return obs_cov, lowerCI, upperCI

		
def bootstrap_mean_95CI(list, reps):
	'''
	calculate the 95% CIs of the mean list of numbers by bootstrapping 'reps' times
	output returned as a tuple: lowerCI, upperCI
	'''
		
	import sys
	sys.path.append( "/Users/Rob/Documents/python_scripts/" )
	from RobPhyloFunctions import bootstrap
	from numpy import mean, median

	observed_mean = mean(list)

	if reps<100:
		print "Less than 100 replications is a bit pointless, so resetting to do 100 replications"
		reps = 100
		
	upperlimit = int(round(reps*0.975))
	lowerlimit = int(round(reps*0.025))
	
	means = []
	
	for i in range(reps):
		means.append(mean(bootstrap(list)))
		
	means.sort()
	
	upperCI = means[upperlimit]
	lowerCI = means[lowerlimit]
	
	return observed_mean, lowerCI, upperCI
	
def bootstrap_median_95CI(list, reps):
	'''
	calculate the 95% CIs of the median list of numbers by bootstrapping 'reps' times
	output returned as a tuple: median, lowerCI, upperCI
	'''
		
	import sys
	sys.path.append( "/Users/Rob/Documents/python_scripts/" )
	from RobPhyloFunctions import bootstrap
	from numpy import mean, median
	
	observed_median = median(list)

	if reps<100:
		print "Less than 100 replications is a bit pointless, so resetting to do 100 replications"
		reps = 100
		
	upperlimit = int(round(reps*0.975))
	lowerlimit = int(round(reps*0.025))
	
	medians = []
	
	for i in range(reps):
		medians.append(median(bootstrap(list)))
		
	medians.sort()
	
	upperCI = medians[upperlimit]
	lowerCI = medians[lowerlimit]
	
	return observed_median, lowerCI, upperCI
	
	
def count_thymine_dimers(sequence):
	'''
	NB sequences must be in UPPERCASE
	take a list of nucleotides, and figure out how many thymine dimers there are
	uses a sliding window of 2 nucleotides
	ignores ambiguity characters, question marks, etc
	strips all gaps (coded as '-') before looking
	
	returns a number (duh)
	'''
	
	seq = sequence.replace('-', '')
	
	pos = 0
	found = 0
	
	for i in range(len(seq)): #longer than I need, but I couldn't figure out how to do it ina while loop
		
		pos = seq.find('TT', pos)
		pos+=1
		
		if pos==0: #it will be zero if the last attempt to find a 'tt' failed
			break
		else:
			found+=1
			
	return found

def count_dimers(sequence, dimers):
	'''
	Input a DNA sequence as a string, and a list of dimers to look for, e.g. ['AA', 'AC', 'AT']
	Output - the number of dimers in the sequence which are found in the list 'dimers'
	overlaps count, e.g. tttt contains 3 tt dimers
	
	case doesn't matter
	
	Gaps are stripped, and should be coded as '-'
	Ambiguity characters are ignored
	'''
		
	#strip gaps
	seq = sequence.replace('-', '')

	#make string uppercase
	seq = seq.upper()
	
	count = 0
	
	#use a sliding window of two, and just count the dimers
	for i in range(len(seq) - 2):
		
		dimer = seq[i:i+2]
		if dimers.count(dimer)>0:
			count+=1
	
	return count

	
def count_pyrimidine_dimers(sequence):
	'''
	Input a DNA sequence as a string
	Output - the number of pyrimidine dimers in the sequence
	Pyrimidine dimers are 'cc' 'ct' 'tc' 'tt'
	overlaps count, e.g. actcta contains 3 dimers
	
	case doesn't matter
	
	Gaps are stripped, and should be coded as '-'
	Ambiguity characters are ignored
	'''
	
	#pyrimidine dimers list
	pd = ['CC', 'CT', 'TC', 'TT']
	
	#strip gaps
	seq = sequence.replace('-', '')

	#make string uppercase
	seq = seq.upper()
	
	count = 0
	
	#use a sliding window of two, and just count the dimers
	for i in range(len(seq) - 2):
		
		dimer = seq[i:i+2]
		if pd.count(dimer)>0:
			count+=1
	
	return count
	
	
	
def calculate_phylo_average(spplist, tree, data):
	'''this function calculates the phylogenetic average of the stuff in data 
	   given MONOPHYLETIC list of species 
	   it does it by successively averaging according to the phylogeny
	   
	   tree should be a tree in biopython format
	   spplist is a list of spp names from the tree as strings
	   data is a dict keyed by species name, with the data as dict entries.
	   
	   '''
	from Bio.Nexus import Nexus, Nodes
	
	#test the monophyly, rootnode will be the common ancestor of all the species in the list
	rootnode = tree.is_monophyletic(spplist)
	if rootnode == -1:
		print "THE INPUT LIST OF TAXA IS NOT MONOPHYLETIC ON THE INPUT TREE"
		print "Input list of taxa:", spplist
	
	else:
					
		#get the list of taxa as a set, in a string
		sppset = str(tree.set_subtree(rootnode))

		#parse the string to remove '[', ']' and 'ImmutableSet' and 'frozenset'
		sppset = sppset.replace('[', '')
		sppset = sppset.replace(']', '')
		sppset = sppset.replace('ImmutableSet', '')
		sppset = sppset.replace('frozenset', '')
		
		#now make it into a sum, by replacing ',' with '+'; '(' with '((', and ')' with ')/2.0)'
		sppsum = sppset.replace(',', '+')
		sppsum = sppsum.replace('(', '((')
		sppsum = sppsum.replace(')', ')/2.0)')
	
	#now get the data for your species list
	#the actual root of all the species is the parent of the current root node
	rootnode = tree.node(rootnode).get_prev()
	for spp in spplist:
		
		#sometimes I use a hack where i put the family name on the front of a species, seperated by '__'. Check for this and remove it if necessary
		if spp.find('__')>0:
			sppstripped = spp.split('__')[1]
		
		else:
			spptripped = spp
		
		sppdata = data[sppstripped]
		
		#now replace that taxon name in the sppsum with the root-to-tip branchlength
		sppsum = sppsum.replace(str(spp), str(sppdata))
		
	#now strip the " ' " from teh sppsum
	sppsum = sppsum.replace("'", '')
	
	#now evaluate that expression you've just made, print it and return it
	sppsum = eval(sppsum)

	print sppset
	print sppsum

	return sppsum
	
def get_sister_pairs(treefilename):
	'''
	This function takes a treefile name (nexus format with taxon list please) as input
	and returns a biopython tree object and a list of comparison objects which look like this
	
	class Comparison:
	def __init__(self, id, spp1name, spp1brlen, spp2name, spp2brlen, brlenratio, sqrtsumbrlen):
		
		self.id 			= id of the root node of the comparison
		self.spp1name 		= spp1name
		self.spp1brlen		= spp1brlen
		self.spp2name 		= spp2name
		self.spp2brlen		= spp2brlen
		self.brlenratio 	= brlenratio (NAN if one of the brlens is zero)
		self.sqrtsumbrlen 	= sqrtsumbrlen (NAN if one of the brlens is zero)

	the return is: tree, list_of_comparisons
		
	'''
	from Bio.Nexus import Nexus, Nodes
	from math import sqrt, log
	import csv
	
	class Comparison:
		def __init__(self, id, spp1name, spp1brlen, spp2name, spp2brlen, brlenratio, sqrtsumbrlen):
			
			self.id 			= id
			self.spp1name 		= spp1name
			self.spp1brlen		= spp1brlen
			self.spp2name 		= spp2name
			self.spp2brlen		= spp2brlen
			self.brlenratio 	= brlenratio
			self.sqrtsumbrlen 	= sqrtsumbrlen
	
	treefile = Nexus.Nexus(treefilename)
	
	tree = treefile.trees[0]
	
	comparisons = []

	#get a list of all the nodes
	all_nodes = tree.all_ids()
	
	#go through each node and get a list of nodes which have sister pairs
	sister_pair_nodes = []
	for node in all_nodes:
		if len(tree.get_taxa(node)) == 2:
			sister_pair_nodes.append(node)
	
	#now for each of those nodes, get the data and put it into a comparison object
	for parent in sister_pair_nodes:
		
		#get the two daughters of each sister pair node
		spp1, spp2 = tree.node(parent).get_succ()
		
		spp1name = tree.get_taxa(spp1)[0]
		spp2name = tree.get_taxa(spp2)[0]
		
		spp1brlen = float(tree.distance(spp1, parent))
		spp2brlen = float(tree.distance(spp2, parent))
		
		#avoid zero branchlengths, just by not including them in the list
		if spp1brlen*spp2brlen>0.0000000:
			
			brlenratio   = log(spp1brlen/spp2brlen)
			sqrtsumbrlen = sqrt(spp1brlen + spp2brlen)

		else:
			brlenratio = 'NAN'
			sqrtsumbrlen = 'NAN'
					
		#put it all in a nice comparison class
		comp = Comparison(parent, spp1name, spp1brlen, spp2name, spp2brlen, brlenratio, sqrtsumbrlen)
		
		comparisons.append(comp)
	
	return tree, comparisons	
	
def garland_test(sqrt_time, trait_differences):
	'''
	Input a pair of lists with corresonding values
	the first contains estimates of the sqrt of time (e.g. sqrt sum of brlens)
	the second contains differences in subsitution rates or traits
	
	This approach uses kendall's tau to test for association between
	
	the sqrt time, and the trait difference.
	
	consider carefully beforehand if and how you want to standardise your trait differences

	it will also save a plot in the current working directory called garland_plot
	
	'''

	from pylab import figure, plot, title, ylabel, xlabel, savefig, rc
	import sys
	sys.path.append( "/Users/Rob/Documents/python_scripts/" )
	import stats



	traits_abs = []
	for thing in trait_differences:
		traits_abs.append(abs(thing))
	
	######################################   GARLAND TEST   ############################################
	#Plot standardised contrast magnitude against the mean sqrt of the sum of the brlens	
	
	tau_garland, garland_p  = stats.kendalltau(sqrt_time, traits_abs)
	
	
	
	print "Garland Test on Contrasts:\n p=%.3f\ntau=%.3f" %(garland_p, tau_garland)

	figure(1)
	rc("font", size = 6)

	plot(sqrt_time, traits_abs, 'ro')
	title('Garland Rates: kendalls_tau=%.3f, p=%.3f' %(tau_garland, garland_p))
	ylabel('absolute standardised contrasts')
	xlabel('sqrt of time')

	figtitle = "garland_plot"
	savefig(figtitle, dpi = 300)


def count_dimers_in_metagenome_data(filename, format, dimers):
	'''
	This function is written to count the number of specified dimers in a given metagenome survey
	Input
		filename: the name of the sequence file with all the dna sequences in
		format:	the format of the sequence file, usually fasta
		dimers: a list of the dimers you're looking for, e.g. if you were searching for thymine dimers it would be ['TT']

	output is prined to the screen, i am too lazy.
	'''

	import sys
	sys.path.append( "/Users/Rob/Documents/python_scripts/" )
	from RobPhyloFunctions import count_dimers
	from Bio import SeqIO
	
	print "\n"
	
	sys.stdout.write("Loading sequences")
	sys.stdout.flush()

	#import the alignemt
	handle = open(filename, "rU")
	
	#make the alignment into a nice list
	seqs = []
	for record in SeqIO.parse(handle, format):
		seqs.append(record.seq.tostring())
			
	#make a dict of the number of thymine dimers
	i = 0
	pyrimidine_dimers = 0
	total_dimers = 0
	total_gc = 0
	total_tc = 0
	total_bases =0
	total_t = 0
	last_percent_done = 0
	for seq in seqs:
	
		percent_done = (100.0*float(i))/float(len(seqs))
	
		if percent_done>last_percent_done:
			sys.stdout.write("\rWorking on Sequences, Percent done: %d\r" %int(percent_done))
			sys.stdout.flush()
			
			
			last_percent_done += 1
		i += 1
		total_dimers += (len(seq) - 1)
		pyrimidine_dimers += count_dimers(seq, dimers)
		a = float(seq.count('A'))
		c = float(seq.count('C'))
		g = float(seq.count('G'))
		t = float(seq.count('T'))
	
		total_bases += a+c+g+t
		total_gc += g+c
		total_tc += t+c
		total_t  += t
	
	prop_dimers = float(pyrimidine_dimers)/float(total_dimers)
	proportion_gc = float(total_gc/total_bases)
	proportion_tc = float(total_tc/total_bases)
	proportion_t  = float(total_t/total_bases)
	print "done calculating stats" 
	
	#print out fam1, fam2, uvdiff, ttdiff, tdiff
	print "\n\n"
	print filename
	print "total number of sequences:", i
	print "total number of bases:", int(total_bases)
	print "proportion of pyrimidine dimers:", prop_dimers
	print "GC%:", proportion_gc
	print "TC%:", proportion_tc
	print "T%", proportion_t
	
	return filename, i, total_bases, prop_dimers, proportion_gc, proportion_tc, proportion_t
	

def get_accession(query, database, rettype):
	#!/usr/bin/env python
	# A short script to download nucleotide sequences from genbank.
	#
	# Copyright (c) 2009, Simon J. Greenhill <simon@simon.net.nz>
	# All rights reserved.
	#
	# Redistribution and use in source and binary forms, with or without
	# modification, are permitted provided that the following conditions
	# are met:
	#
	# 1. Redistributions of source code must retain the above copyright notice,
	#    this list of conditions and the following disclaimer.
	#
	# 2. Redistributions in binary form must reproduce the above copyright
	#    notice, this list of conditions and the following disclaimer in the
	#    documentation and/or other materials provided with the distribution.
	#
	# 3. Neither the name of genbank-download nor the names of its contributors
	#    may be used to endorse or promote products derived from this software
	#    without specific prior written permission.
	#
	# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
	# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
	# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
	# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
	# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
	# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
	# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
	# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
	# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
	# POSSIBILITY OF SUCH DAMAGE.
	"""
	Returns a nucleotide sequence from genbank.
	
	:param query: the accession number
	:param database: the database to use (default=nucleotide)
	:param rettype: the return type for the sequence (native,fasta,gb,xml)
	
	:return: text of sequence in requested `rettype`
	:rtype: string
	
	"""
	
	import urllib
	
	_toolname = 'get_accession'
	_email = "rob.lanfear@gmail.com"
	
	params = {
		'db': database,
		'tool': _toolname,
		'email': _email,
		'id': query,
		'rettype': rettype,
	}
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
	url = url + urllib.urlencode(params)
	data = urllib.urlopen(url).read()
	return data


def get_accession_length(accession):
	"""	get the length of an accession from genbank"""
	import sys
	sys.path.append( "/Users/Rob/Documents/python_scripts/" )
	from RobPhyloFunctions import get_summary
	import re
	acc = get_summary(accession, 'nucleotide', 'xml')
	#length looks like this: <Item Name="Length" Type="Integer">1590</Item>
	length = re.search('<Item Name="Length" Type="Integer">(\S+)</Item>', acc)
	return length.group(1)
	

def get_taxid(term):
	#get a taxid from entrez taxonomy, by searching with a name
	import urllib
	import re
	_toolname = 'get_taxid'
	_email = "rob.lanfear@gmail.com"
	params = {
		'db': 'taxonomy',
		'tool': _toolname,
		'email': _email,
		'term': term,
		'rettype': 'xml',
	}
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
	url = url + urllib.urlencode(params)
	data = urllib.urlopen(url).read()
	taxid = re.search('<Id>(\S+)</Id>', data).group(1)
	return taxid


def get_summary(query, database, rettype):
	#this is copied from Simon Greenhill's get_accession function, see get_accession for copyright, above
	import urllib
	_toolname = 'get_summary'
	_email = "rob.lanfear@gmail.com"
	params = {
		'db': database,
		'tool': _toolname,
		'email': _email,
		'id': query,
		'rettype': rettype,
	}
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'
	url = url + urllib.urlencode(params)
	data = urllib.urlopen(url).read()
	return data

def get_all_GIs(taxid):
	"""
	Get all GIs for a taxon ID, that match that ID
	:param taxid: a genbank taxnomy ID as a string
	"""
	import urllib, re
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?&retmax=200&term=txid%s[Organism:subtree]&' %(taxid)
	params = {
		'db': 'nucleotide',
		'rettype': 'xml',
	}
	url = url + urllib.urlencode(params)
	data = urllib.urlopen(url).read()
	GIs = re.findall("<Id>(\S+)</Id>", data)
	return GIs

def search_all_GIs(taxid, searchterm):
	"""
	Get all GIs for a taxon ID, that match that ID, and the search term
	:param taxid: a genbank taxnomy ID as a string
	:param searchterm: a string, of 1 or more words
	N.B. Don't use spaces in the search term, use commas instead.	
	"""
	import urllib, re
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?&retmax=200&term=txid%s[Organism:subtree]%s&' %(taxid, searchterm)
	params = {
		'db': 'nucleotide',
		'rettype': 'xml',
	}
	url = url + urllib.urlencode(params)
	data = urllib.urlopen(url).read()
	GIs = re.findall("<Id>(\S+)</Id>", data)
	return GIs


def get_longest_accession2(taxid, names_list):
	'''
	like the other get_longest_accession, but a bit more thorough, and uses the proper channels at genbank
	
	searches entrez nucleotide (note, NOT entrez genome) for genes
	
	input:
		taxid: string, corresponds to genbank taxid
		names_list: list of strings, corresponds to possible names of genes
		
	output:
		GI corresponding to longest gene matching input criteria
		if >1 gene are equally long, then a random one is picked
	'''

	import sys
	sys.path.append( "/Users/Rob/Documents/python_scripts/" )
	from RobPhyloFunctions import get_summary, get_all_GIs, search_all_GIs
	import re
	import random
	from xml.dom import minidom

	#1. Get all GIs for that taxid
	GIs = []
	for name in names_list:
		name = name.replace(" ", ",") #get rid of commas, since efetch doesn't like em. 
		temp = search_all_GIs(taxid, name)
		for thing in temp:
			GIs.append(str(thing))
	
	#2. Get the summaries of all the GIs
	GIdict = {}
	
	summary = get_summary(','.join(GIs), 'nucleotide', 'xml') #send in the big list of all GIs
	
	#now parse the GIs
	tempf = open("temp.xml", 'w')
	tempf.write(summary)
	tempf.close()
	
	xmldoc = minidom.parse("temp.xml")
	
	start = xmldoc.childNodes[1]
	nodelist = start.childNodes
	for thing in nodelist:
		if str(thing).count("DocSum")>0:
			GIsummary = thing.toxml()
			GI = str(int(re.search('<Item Name="Gi" Type="Integer">(\S+)</Item>', GIsummary).group(1)))
			GIdict[GI] = GIsummary
	
	#3. Throw out GIs whose name field doesn't match what you're looking for
	goodGIdict = {}
	for GI in GIdict:
		GIname = re.search('<Item Name="Title" Type="String">(.*)</Item>', GIdict[GI]).group(1)
		total_count = 0
		for thing in names_list:
			total_count += GIname.count(thing)
		if total_count > 0: #if you didn't find any matches with the name you wanted
			goodGIdict[GI] = GIdict[GI]

	#4. pick the longest accession that's left
	longest = -999
	accessions = []
	lengthdict = {}

	if goodGIdict: #if there's anything left
		for GI in goodGIdict:
			length = int(str(re.search('<Item Name="Length" Type="Integer">(\S+)</Item>', goodGIdict[GI]).group(1)))
			lengthdict[length] = GI
		lengths = lengthdict.keys()
		lengths.sort()
		lengths.reverse()
		longest = lengths[0]
		longest_acc = lengthdict[longest]		

	else:
		longest_acc = 'x'
		
	#if you found>1 accession, randomly choose one
	if accessions:
		longest_acc = random.choice(accessions)
	
	return longest_acc
	
	
			


def get_longest_accession(taxid, names_list, genome):
	'''
	This function is designed for trawling genbank to build supermatrices 
	it takes as input a taxid from genbank taxonomy, and a list of names to look for, and which genome the gene is on
	(e.g. the list of names might look like this CYTB_list = ['cytb', 'CYTB', 'cytochrome b', 'CYTOCHROME B', etc])
	and a genome (i.e. "mitochondrial", "nuclear", or "chloroplast")
	the function also has an inbuilt dict of these lists, for which you can just give the name
	
	it returns the longest accesssion it can find, for which the name of the accession matches
	one of the terms in the list you have supplied
	
	if it can't find that accession, it returns 'x'
	'''
	
	import urllib2
	import re
	import random
	random.seed = 127 #for repeatability...
	
	if names_list == "CYTB":
		names_list = ["mitochondrial genome", "mitochondrial complete genome", "Mitochondrial Genome", "Mitochondrial genome", "CYTB", "cytb", "cytochrome b", "Cytochrome b", "Cytochrome B", "cytochrome oxidase b", "Cytochrome Oxidase B"]
	if names_list == "COI":
		names_list = ["mitochondrial genome", "mitochondrial complete genome", "Mitochondrial Genome", "Mitochondrial genome", " COI ", " CO1 ", "(COI)", "(CO1)", "cytochrome oxidase 1 ", "Cytochrome Oxidase 1 ", "cytochrome oxidase subunit 1 ", "Cytochrome Oxidase Subunit 1 ", "Cytochrome Oxidase I ", "cytochrome oxidase subunit I ", "Cytochrome Oxidase Subunit I ", "cytochrome oxidase subunit I "]
	if names_list == "12S":
		names_list = ["mitochondrial genome", "mitochondrial complete genome", "Mitochondrial Genome", "Mitochondrial genome", "12S", "large subunit ribosomal RNA"]
	if names_list == "ND2":
		names_list = ["mitochondrial genome", "mitochondrial complete genome", "Mitochondrial Genome", "Mitochondrial genome", "ND2", "NADH2", "NADH dehydrogenase subunit 2"]
	if names_list == "CMOS":
		names_list = ["CMOS", "oocyte maturation factor", "c-mos", "C-MOS", "C-Mos"]
	if names_list == "RAG1":
		names_list = ["RAG1", "RAG-1", "recombination activating protein 1"]
	if names_list == "MYO":
		names_list = ["(myo)", "(MYO)", " myo ", " MYO ", "myoglobin gene", " Myo ", "(Myo)"] 
		
	#build the url for entrez nuc
	urlend = '%5BOrganism%3Aexp%5D&go=Go&dispmax=1000000'
	urlstart = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=nuccore&cmd=&term=txid%s' %(taxid)
	url = ''.join([urlstart, urlend])


	###################################################################
	#if it's a mito genome sequence, look at the genome database first
	accession = 'x'
	if genome=="mitochondrial":

		#build the URL for entrez genome
		urlend_g = "%5BOrganism%3Aexp%5D&go=Go&dispmax=1000000"
		urlstart_g = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=genome&cmd=Search&dopt=DocSum&term=txid%s" %(taxid)
		url_g = ''.join([urlstart_g, urlend_g])

		#get the url
		entrez_g_result = urllib2.urlopen(url_g).read()
	
		#save to a temporary file
		tempf = open("tempg.xml", "w")
		tempf.write(entrez_g_result)
		tempf.close()
		
		#read in tempf as lines
		file = open("tempg.xml", "r")
		lines = file.readlines()
			
		accessions = [] #an empty list of accessions
		for line in lines:
			#the accession looks like this: >NC_011493</a></td><td align
			a = re.search(">NC_(\d+)</a></td><td ", line)
			if a:
				accessions.append(''.join(["NC_", a.group(1)]))
	
		#randomly choose one of the whole-genome accessions
		if accessions:
			accession = random.choice(accessions)

	###################################################################
	#if you didn't get an accession from the genome database, look in the nuccore	
	if accession == 'x':
		#get the url
		entrez_result = urllib2.urlopen(url).read()
	
		#save to a temporary file
		tempf = open("temp.xml", "w")
		tempf.write(entrez_result)
		tempf.close()
		
		#read in tempf as lines
		file = open("temp.xml", "r")
		lines = file.readlines()
		
		#entrez results start with "<div class="docsum" "
		results = []
		for line in lines:
			if line.startswith('<div class="docsum"'):
				#see if that entry contains any of the names of interest			
				total_count = 0
				for thing in names_list:
					total_count += line.count(thing)
				if total_count>0:
					results.append(line)
		
		#parse results
		longest = -999
		accessions = []
		for line in results:
			#get the length. lengths look like this: >1,041 bp
			l = re.search('>(\S+) bp', line)
			length = int(l.group(1).replace(",", ""))
			#the accession looks like this: href="/nuccore/EU236668.1" ref=
			a = re.search('href="/nuccore/(\S+)" ref=', line)
			#if it's equal record the id and add to the list
			if length==longest:
				accessions.append(a.group(1))
			#if it's longer, re-set the list and the longest
			if length>longest:		
				#re-set longest to be length
				longest = length
				accessions = []
				accessions.append(a.group(1))
	
		
		if accessions:
			accession = random.choice(accessions)

	return accession	

def get_taxonomy(accession):
	'''
	This function takes as input a genbank accession (for the nucleotide database only)
	it returns three things: taxomony(dictionary), binomial(string), TaxonID(string)
	the taxonomy dictionary looks like this
		taxonomy{
			phylum: ??
			subphylum: ??
			superclass: ??
			class: ??
			superorder: ??
			order: ??
			suborder: ??
			family: ??
			subfamily: ??
			genus: ??
			species: ??
		}
	if there's data for that taxonomic level in genbank, then the appropriate level will be filled in,
	otherwise, they're left as question marks "??".
	'''
	
	import urllib2
	import re
	
	#taxonomy is a dict with levels (e.g. phylum) as the keys, and names as entries (e.g. chordata)
	#initialised with empty question marks
	taxonomy = Taxonomy = {'phylum': '??', 'subphylum': '??', 'superclass': '??', 'class': '??', 'superorder': '??', 'order': '??', 'suborder': '??', 'family': '??', 'subfamily': '??', 'genus': '??'}
	

	#first go look up the accession on entrez nucleotide (assuming of course that it's a nucleotide accession...)
	url = 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nuccore&id=%s' %(accession)
	entrez_nuc = urllib2.urlopen(url).read() #this is the dirty looking xml file from entrez
	
	#extract the TaxonID from the horrible looking xml file 
	m = re.search('/Taxonomy/Browser/wwwtax.cgi\?id=(\d+)', entrez_nuc)
	TaxID = m.group(1)
	
	#now use the taxon ID to go and get the full lineage from entrez taxonomy
	url = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=%s&lvl=3&keep=0&srchmode=1&unlock&lin=f' %(TaxID)
	entrez_tax = urllib2.urlopen(url).read() #this is the dirty looking xml file from entrez

	#print "length of entrez_tax =", len(entrez_tax)
	
	#now extract all the lineage data from entrez_tax
	for level in taxonomy:
		#this gets the different levels from the entrez_tax html
		m = re.search('"%s">(\S+)</a>' %(level), entrez_tax)
		#if the search stuff worked, and you found that level...
		if m:
		#only do this bit if you actually found that part of the taxonomy, otherwise just leave it blank.
		#this prevents the script unnecessarily hanging up.
			taxonomy[level] = m.group(1)

	#now get the binomial (or in some cases its a trinomial, but whatever)
	binomial = re.search(r">Taxonomy browser ((.+))<", entrez_tax)
	binomial = binomial.group(1)[1:-1]
	species = binomial.split()[1]
	
	#put the species name in the taxonomy (the reason I add it in here, is it's not listed in the entrez_tax listing I used above
	taxonomy['species'] = species
	
	#print taxonomy
	return taxonomy, binomial, TaxID

def find_and_replace(inputfileName, outputfileName, sourceText, replaceText):
	'''
	NB. This is a slightly modified version of a find-replace function I found here:
	http://www.linuxquestions.org/questions/programming-9/python-find-defined-text-string-in-a-file-and-replace-the-whole-line-186924/
	I have not been able to find out who wrote the original, so I don't know how to cite this properly.
	If it was you, please tell me, and I'll happily cite you appropriately.
	
	If inputfileName and outputfileName are different, you get a new file
	if they're the same, you overwrite the old input file.
	'''
	file = open(inputfileName, "r")
	text = file.read()
	file.close()
	file = open(outputfileName, "w")
	file.write(text.replace(sourceText, replaceText)) #replaces all instances of our keyword
	file.close()
 
def add_systematic_rate_variation(brlen_tree, labelled_tree, rates_dict):
	'''Takes a tree with brlens, and one with labels, as inputs
	as well as dictionary which relates labelled to the rate multiplier you want
	Then it goes through and multiples the labelled branches by
	the appropriate rate multiplier
	NB: the labelled and brlen trees should be identical in structure
	'''
	######################################################
	#make a dict of branch labels, keyed by branch_id
	branch_labels = {}
	is_labelled = 0
	branchnumber = 1
	i = 0
	for i in range(len(labelled_tree)):
		thing = labelled_tree[i]
		#reset the branch, just in case
		branch_id = ''
		if thing == '#':
			#flag the branch labels as they appear, this flag
			#gets reset to zero after you find the next branch
			is_labelled = 1
			#extract the label - this will ONLY work for single-digit labels
			branch_label = '#%s' %(labelled_tree[i+1])
		if thing == ',' or thing == ')':
			#it's a branch
			#the list of branch labels is added to the new tree, just to demarcate which branches are which
			branch_id = "%s%d" %('##branch_', branchnumber)
			if is_labelled == 1:
				branch_labels[branch_id] = branch_label
			is_labelled = 0
			branchnumber += 1
	######################################################
	#make a dict of branch lengths, keyed by branch_id
	branch_lengths = {}
	is_brlen = 0
	branchnumber = 1
	i = 0
	for i in range(len(brlen_tree)):
		thing = brlen_tree[i]
		#reset the branch, just in case
		branch_id = ''
		if thing == ':':
			#flag the branch labels as they appear, this flag
			#gets reset to zero after you find the next branch
			is_brlen = 1
			#extract the label - this will ONLY work for single-digit labels
			brlen = []
			j = i + 1
			while (brlen_tree[j].isdigit() or brlen_tree[j] == '.' or brlen_tree[j]==' '):
				brlen.append(brlen_tree[j])
				j = j+1
			brlen = float(''.join(brlen))
		if thing == ',' or thing == ')':
			#it's a branch
			#the list of branch labels is added to the new tree, just to demarcate which branches are which
			branch_id = "%s%d" %('##branch_', branchnumber)
			branch_lengths[branch_id] = brlen
			branchnumber += 1	
	######################################################
	#multiply the brlens by the appropriate value in the rate dict
	for branch in branch_lengths:
		if branch_labels.has_key(branch):
			branch_length = branch_lengths[branch]
			label = branch_labels[branch]
			branch_lengths[branch] = branch_lengths[branch]*rates_dict[label]
	######################################################
	#use the labelled tree to make a new tree, with labels and brlens
	#first make a new tree as a list, and a list of branches
	#in the new tree, the each branch will be labelled like this '##branch_x'
	#where x is an integer. Hopefully this won't conflict with any taxon labels...
	newtree = []
	branches = []
	# a dict of branch labels, keyed by branch_id
	is_labelled = 0
	branchnumber = 1
	i = 0
	for i in range(len(labelled_tree)):
		thing = labelled_tree[i]
		#reset the branch, just in case
		branch_id = ''
		if thing == '#':
			#flag the branch labels as they appear, this flag
			#gets reset to zero after you find the next branch
			is_labelled = 1
		if thing == ',' or thing == ')':
			#it's a branch
			#the list of branch labels is added to the new tree, just to demarcate which branches are which
			branch_id = "%s%d" %('##branch_', branchnumber)
			branches.append(branch_id)
			if is_labelled == 1:
				branch_labels[branchnumber] = branch_label
			branchnumber += 1
			is_labelled = 0
		newtree.append(branch_id)
		if is_labelled == 0:
			#in this way, the new tree is free of all labels
			newtree.append(thing)
	#################################################################
	#now make a new tree, with labels and branch lengths
	#now replace all branches in the tree with the appropriate label
	#if the branch isn't in branchdict, that means it doesn't get a lable (and PAML will assign it to rate category zero)
	for branch in branches:
		branch_length = str(branch_lengths[branch])
		if branch_labels.has_key(thing):
			newlabel = ''.join([branch_labels[thing], ' :', branch_length])
		else:
			newlabel = ' :%s' %(branch_length)
		newtree[newtree.index(branch)] = newlabel
	#now turn the newtree list into a phylip looking tree
	outputtree = ''.join(newtree)

	return outputtree
	
def extract_likelihoods_from_PAML_output(file):
	'''
	Input a PAML file,
	get a list of all the likelihoods in that file
	'''
	
	lines = file.readlines()
	likelihoods = []
	for line in lines:
		if line.startswith("lnL(ntime:"):
			temp = line.split("):")[1]
			temp = temp.strip()
			likelihood = float(temp.split(' ')[0])
			likelihoods.append(likelihood)
	return likelihoods
	
	
def extract_relrates_from_PAML_output(file, N):
	'''
	Input a PAML file,
	get a list of all the relative rates in that file
	you need to specify the filename, and the number of rate parameters you're expecting
	N.B. This function is only appropriate for output files with 1 tree analysed
	'''
	lines = file.readlines()
	rates = []
	for line in lines:
			if line.startswith("rates for branches:"):
				temp = line.split()
				for i in range(N):
					rates.append(float(temp[3+i]))
	return rates
	