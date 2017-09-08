#!/usr/bin/env python

from math import exp,log,sqrt
import numpy as np
import dendropy

#assumes dt sufficiently small so no more than 1 event will happen
def calcBirthProb(brate,dt):
	return 1.0-exp(-brate*dt)

def calcDeathProb(drate,dt):
	return 1.0-exp(-drate*dt)

	# for coalescing lineages

def mergeLineages(lin1,lin2):
	node=dendropy.Node()
	node.gens=0
	node.edge_length=0.0
	node.add_child(lin1)
	node.add_child(lin2)
	return node

def mergeLineagesMT(lin1,lin2):
	node=mergeLineages(lin1,lin2)
	node.itype=lin1.itype
	node.muts=[]
	assert lin1.itype==lin2.itype
	return node

def extendLineage(lin,dt,gens):
	lin.edge_length=lin.edge_length+dt
	lin.gens=lin.gens+gens
	return lin 

def initLineage():
	newnode=dendropy.Node()
	newnode.edge_length=0.0
	newnode.gens=0
	return newnode

def initLineageMT(itype):
	newnode=initLineage()
	newnode.itype=itype # type at the top, not bottom!
	newnode.muts=[] # mutations on the branch
	return newnode

def addMutationToLineage(lineage,prevtype):
	nexttype=lineage.itype
	lineage.muts.append((nexttype,lineage.edge_length))
	lineage.itype=prevtype


def printLineages(lineages):
	for lineage in lineages:
		t=dendropy.Tree(seed_node=lineage)
		print str(t)+';\n'


def scaleByGenerations(tree):

	gentree=dendropy.Tree(tree)
	for node in gentree.postorder_node_iter():
		if node.edge_length!=None:
			node.edge_length=node.gens
	return gentree

def sampleMultiVariateHyperGeometric(popvec,nSamples):
	totalleft=sum(popvec)
	sampleleft=nSamples
	assert nSamples<=totalleft

	sample=[]
	for i in range(len(popvec)-1):
		totalleft-=popvec[i]
		if sampleleft==0:
			sample.append(0)
		else:
			smpl=np.random.hypergeometric(popvec[i],totalleft,sampleleft)
			sample.append(smpl)
			sampleleft-=smpl
	#last category
	sample.append(sampleleft)
	return sample

#converts a list of counts into a longer list where each index is repeated the number of times it occurs in the list of counts
#[2,1,3] -> [0,0,1,2,2,2]
def serializeCounts(counts):
	outlist=[]
	for i in range(len(counts)):
		for j in range(counts[i]):
			outlist.append(i)
	return outlist


class Deme:

	weight=1.0 # relative weight compared to sister demes - in a population tree, the ancestral pop. gets divided according to those weights
	Nstart=100 # number at start
	Nend=None #number after simulating forward
	Nsampled=0 # number of those that we sample at the end
	tStart=0
	tEnd=10.0
	dt=1.0
	brate=1.0
	drate=0.0

	pbirth=None
	pdeath=None

	simSteps=[]
	popSizeSteps=[] #population at every step= simSteps[cur][1]+2*simSteps[cur][2]
	bottomLineages=[]



	def __init__(self,brate,drate,dt,ts,te):
		self.brate=brate
		self.drate=drate
		self.tStart=ts
		self.tEnd=te
		self.dt=dt
		self.pbirth=calcBirthProb(brate,dt)
		self.pdeath=calcDeathProb(drate,dt)
		self.simSteps=[]
		self.popSizeSteps=[]
		self.bottomLineages=[]

	def resetSimulation(self):
		simSteps=[]


	def simulateForward(self,N):
		self.Nstart=N
		Ncur=N
		nsteps=int((self.tEnd-self.tStart)/self.dt)

		#print 'tstart='+str(self.tStart)+'tend='+str(self.tEnd)+'nsteps='+str(nsteps)

		for i in range(nsteps):
			Ncur=self.simStepForward(Ncur)
			self.popSizeSteps.append(Ncur)

		self.Nend=Ncur
		#print self.Nend




	def simStepForward(self,N):
		distr=[self.pdeath,1.0-self.pbirth-self.pdeath,self.pbirth]
		events=list(np.random.multinomial(N,distr))
		self.simSteps.append(events)
		Nnew=events[1]+2*events[2]
		return Nnew

	def simStepBackward(self,step,lineages):
		#print "base simStepBackward"

		curevents=self.simSteps[step]

		njustborn=curevents[2]*2
		justbornsampled={}
		#print self.popSizeSteps[step]
		#print len(lineages)
		ancestors=random.sample(xrange(self.popSizeSteps[step]),len(lineages))

		lincur=[]

		for i in range(len(lineages)):
			anc=ancestors[i]
			if anc<njustborn:
				extendLineage(lineages[i],self.dt,1)# add one generation to cell generation branch length
				justbornsampled[anc]=lineages[i]
				sisterid=(anc/2)*2+(anc+1)%2
				if sisterid in justbornsampled:
					newlin=mergeLineages(lineages[i],justbornsampled[sisterid])
					lincur.append(newlin)
					del justbornsampled[sisterid]
					del justbornsampled[anc]
			else:
				extendLineage(lineages[i],self.dt,0)
				lincur.append(lineages[i])

		for i in justbornsampled:
			lincur.append(justbornsampled[i])

		return lincur



	def simulateBackward(self,lineages):

		nsteps=len(self.simSteps)

		#print 'BACKWARD: tstart='+str(self.tStart)+'tend='+str(self.tEnd)+'nsteps='+str(nsteps)
		lincur=lineages
		for i in range(nsteps-1,-1,-1):
			lincur=self.simStepBackward(i,lincur)

		return lincur

	def sampleAtBottom(self,nSamples):
		self.Nsampled=nSamples
		lineages=[]
		for i in range(nSamples):
			newnode=initLineage()
			lineages.append(newnode)
		self.bottomLineages=lineages
		return lineages


