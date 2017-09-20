#!/usr/bin/env python

import dendropy
import sys
from math import exp,log,sqrt
import numpy as np
import scipy
import random
import copy
import pdb


from class_Deme_jakub import *
from class_MultiTypeDeme_jakub import *
from class_History_jakub import *


def annotate_tree_with_state_changes(atree):

	tree=dendropy.Tree(atree)

	for node in tree.postorder_node_iter():
		if len(node.muts)>0:
			statestring=str(node.itype)
			for (state,time) in node.muts:
				statestring+=str(state)
			for (state,time) in node.muts:
				timestring+=','
				timestring+=str(time)
			node.annotations.add_new('types',statestring)
			node.annotations.add_new('change_times',timestring[1:])

	root=tree.seed_node
	root.annotations.add_new('start_type',str(root.itype))

	return tree



def sample_mut_tree(atree,mut_rate,add_driver_muts=True):

	tree=dendropy.Tree(atree)

	for node in tree.postorder_node_iter():
		rate=node.edge_length*mut_rate
		nmuts=np.random.poisson(rate)
		if add_driver_muts:
			nmuts+=len(node.muts)
		node.edge_length=nmuts

	return tree


def sample_muts_on_lineages(lineages,mut_rate,add_driver_muts=True):
	linsimple=convert_to_simple_lineage_list(lineages)
	linsastrees=list(map(lambda x:dendropy.Tree(seed_node=x),linsimple))
	linmuts=list(map(lambda x: sample_mut_tree(x,mut_rate),linsastrees))
	return linmuts




#converts a tree into a site-frequency spectrum. Assumes that the branch length is proportional to the number of mutations on the branch
def tree_to_SFS(tree,normalize=True):
	nleaves=len(tree.leaf_nodes())
	sfs=[0.0 for i in range(nleaves+1)]

	for node in tree.postorder_node_iter():
		if node.is_leaf():
			node.ndesc=1
		else:
			node.ndesc=sum([nd.ndesc for nd in node.child_nodes()])

		sfs[node.ndesc]+=node.edge_length

	if normalize:
		total=sum(sfs)
		sfsout=[freq/total for freq in sfs]
		return sfsout

	return sfs

#convert to a one-level list of lineages
def convert_to_simple_lineage_list(lineages):
	alllineages=[]
	if isinstance(lineages[0],list):# it's a list of lists (each corresponds to starting state) so we will convert it to a simple list
		for linlst in lineages:
			for lin in linlst:
				alllineages.append(lin)
	else:
		alllineages=lineages
	return alllineages

def nodes_to_trees(lineages):
	alllineages=convert_to_simple_lineage_list(lineages)
	trees=list(map(lambda x: dendropy.Tree(seed_node=x),alllineages))
	return trees

def trees_to_nodes(trees):
	lins=list(map(lambda x: x.seed_node,trees))
	return lins

def lineages_to_SFS(lineages,nsampled=None,normalize=True):
	alllineages=convert_to_simple_lineage_list(lineages)
	#check the total num of sampled individuals if not supplied
	if nsampled==None:
		nsampled=0
		for lin in alllineages:
			tree=dendropy.Tree(seed_node=lin)
			nleaves=len(tree.leaf_nodes())
			nsampled+=nleaves

	# main functionality here
	sfs=[0.0 for i in range(nsampled+1)]
	print("lensfs="+str(len(sfs)))

	for lineage in alllineages:
		tree=dendropy.Tree(seed_node=lineage)
		sfslin=tree_to_SFS(tree,normalize=False)
		for i in range(len(sfslin)):
			sfs[i]+=sfslin[i]

	if normalize:
		total=sum(sfs)
		sfsout=[freq/total for freq in sfs]
		return sfsout

	return sfs



def print_lineages(lineages):
	alllineages=convert_to_simple_lineage_list(lineages)
	genenum=1
	outstr=''
	for lineage in alllineages:
		if(isinstance(lineage,dendropy.Node)):
			tree=dendropy.Tree(seed_node=lineage)
		else:
			tree=lineage
		tree.randomly_assign_taxa()
		for leaf in tree.leaf_nodes():
			leaf.taxon.label='T'+str(genenum)
			genenum+=1
		outstr+=(tree.as_string(schema='newick',
			suppress_annotations=False)+';\n')

	return outstr

# Samples a genealogy from the subsampled process.
# INPUT: Simulation parameters:
# brate - a list of birth rates for each state
# drate - a list of death rates for each state
# dt - time step
# ts - start time (default=0.)
# te - end time
# transmat - transition matrix for per-time-unit, per-cell transition probabilities between states 
# startdist - a list with the number of cells in each state at the start of the simulation
# n_samples_bottom - the number of individuals sampled at end time (if -1, sample all of the population - won't work for large populations)
# mut_rate - the expected number of neutral mutations per time unit (the number of mutations is Poisson-distributed) 

def sample_genealogy(brate,drate,dt,ts,te,transmat,startdist,n_samples_bottom,mut_rate,normalize=True):
	deme=MultiTypeDeme(brate,drate,dt,ts,te,transmat)#set up the class params

	deme.simulateForward(startdist)
	if n_samples_bottom==-1:
		n_samples_bottom=sum(deme.tdEnd) #if not supplied, sample all of the population

	sampledlineages=deme.sampleAtBottom(n_samples_bottom) # set up subsampled lineages from the bottom to sample upwards

	lins=deme.simulateBackward(sampledlineages) 

	treeswithmuts=sample_muts_on_lineages(lins,mut_rate) # place neutral muts on branches (Poisson-distributed)
	linswithmuts=trees_to_nodes(treeswithmuts)
	sfs=lineages_to_SFS(linswithmuts,None,normalize)
	return linswithmuts,sfs




#--------------Test functions-----------

def sample_genealogy_test():
	lins,sfs=sample_genealogy(brate=[0.3,0.4],drate=[0.29,0.29],dt=0.1,ts=0.0,
		te=50,transmat=[[0.999,0.001],[0.0,1.0]],startdist=[100,0],n_samples_bottom=20,mut_rate=0.1,normalize=False)
	print(print_lineages(lins))
	print(sfs)

def multideme_test():
	deme=MultiTypeDeme(brate=[0.3,0.4],drate=[0.29,0.29],dt=0.1,ts=0.0,
		te=50,transmat=[[0.999,0.001],[0.0,1.0]])
	deme.simulateForward([100,0])
	sampledlineages=deme.sampleAtBottom(20)
	lins=deme.simulateBackward(sampledlineages)

	#tree=dendropy.Tree(seed_node=lins[0][0])
	
	#tree.randomly_assign_taxa()
	
	#mtree=sample_mut_tree(tree,0.1)
	#print "muttree="+str(mtree)
	#print "SFS="+str(tree_to_SFS(tree))
	#print str(tree)+';'

	print "SFS2="+str(lineages_to_SFS(lins,normalize=False))

	
	#print print_lineages(lins)
	#linsflat=convert_to_simple_lineage_list(lins)
	#ltrees=nodes_to_trees(linsflat)
	#ltreesproc=map(annotate_tree_with_state_changes,ltrees)

	#print print_lineages(ltreesproc)



def tree_test_simple():
	root=dendropy.Node()
	child=dendropy.Node()
	root.add_child(child)
	root.deme=Deme(3.0,0.0,0.1,0.0,5.0)
	child.deme=Deme(1.0,0.9,0.1,5.0,6.0)
	child.deme.Nsampled=10000
	tree=dendropy.Tree(seed_node=root)
	hist=History(tree)
	lins=hist.simulate(1)
	
	assert len(lins)==1
	tree=dendropy.Tree(seed_node=lins[0])
	gentree=scaleByGenerations(tree)
	tree.randomly_assign_taxa()
	#print str(tree)+';'
	gentree.randomly_assign_taxa()
	print str(gentree)+';'

def tree_test_bifurcation():
	root=dendropy.Node()
	child=dendropy.Node()
	root.add_child(child)
	child2=dendropy.Node()
	root.add_child(child2)

	root.deme=Deme(3.0,0.0,0.1,0.0,5.0)
	child.deme=Deme(1.0,0.9,0.1,5.0,6.0)
	child2.deme=Deme(1.0,0.9,0.1,5.0,6.0)
	child.deme.Nsampled=1000
	child2.deme.Nsampled=1000


	tree=dendropy.Tree(seed_node=root)

	hist=History(tree)
	lins=hist.simulate(1)
	
	assert len(lins)==1
	tree=dendropy.Tree(seed_node=lins[0])
	gentree=scaleByGenerations(tree)
	tree.randomly_assign_taxa()
	print str(tree)+';'
	gentree.randomly_assign_taxa()
	print str(gentree)+';'



#deme_test()
#tree_test_bifurcation()
#multitype_deme_test()
#multideme_test()
sample_genealogy_test()



