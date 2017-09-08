#!/usr/bin/env python

import copy
import random

from class_Deme_jakub import * 

class MultiTypeDeme(Deme):

	tdStart=[100,0]
	tdCur=[]
	brate=[1.0,1.2]
	drate=[0.0,0.0]
	nStates=2
	TransMatrix=[[0.9,0.1],[0.0,1.0]]

	def __init__(self,brate,drate,dt,ts,te,transmat):
		assert len(brate)==len(drate)
		self.nStates=len(brate)
		self.brate=list(brate)
		self.drate=list(drate)
		self.tStart=ts
		self.tEnd=te
		self.dt=dt
		self.pbirth=map(lambda x: calcBirthProb(x,dt),brate)
		self.pdeath=map(lambda x: calcDeathProb(x,dt),drate)
		self.simBDSteps=[]
		self.simTrSteps=[]
		self.popSizeSteps=[]
		self.bottomLineages=[]
		self.TransMatrix=copy.deepcopy(transmat)


	# does the birth-death part of simStepForward - analogous to the original
	# typedist= numbers for each state
	def simBDStepForward(self,typedist):
		assert len(self.pdeath)==len(self.pbirth)
		distr=[[self.pdeath[i], 1.0-self.pbirth[i]-self.pdeath[i], self.pbirth[i]] for i in range(self.nStates)]
		events=[ list( np.random.multinomial(typedist[i],distr[i])) for i in range(self.nStates)]
		self.simBDSteps.append(events)
		typedistNew=[events[i][1]+2*events[i][2] for i in range(self.nStates)]
		return typedistNew

	def simTrStepForward(self,typedist):

		events=[np.random.multinomial(typedist[i],self.TransMatrix[i]) for i in range(self.nStates)]
		self.simTrSteps.append(list(events))
		typedistNew=sum(events)
		return typedistNew
		

	def simTrStepBackward(self,step,lineages):
		curevents=self.simTrSteps[step]
		newlineages=[ [] for i in range(self.nStates)]
		for i in range(self.nStates):
			possibleancs=[ curevents[k][i] for k in range(self.nStates) ]
			anctypes=sampleMultiVariateHyperGeometric(possibleancs,len(lineages[i]))
			ancestors=serializeCounts(anctypes)
			random.shuffle(ancestors)
			#ancestors now has the ancestral state of every lineage in lineages[i]
			for j in range(len(lineages[i])):
				if ancestors[j]!=i:
					addMutationToLineage(lineages[i][j],ancestors[j])
				newlineages[ancestors[j]].append(lineages[i][j])

		return newlineages

	def simBDStepBackward(self,step,lineages):
		
		curevents=self.simBDSteps[step]
		lincur=[ [] for i in range(self.nStates)]
		
		for i in range(self.nStates):

			njustborn=curevents[i][2]*2
			justbornsampled={}
		
			ancestors=random.sample(xrange(curevents[i][1]+njustborn),len(lineages[i]))

		

			for j in range(len(lineages[i])):
				anc=ancestors[j]
				if anc<njustborn:
					extendLineage(lineages[i][j],self.dt,1)# add one generation to cell generation branch length
					justbornsampled[anc]=lineages[i][j]
					sisterid=(anc/2)*2+(anc+1)%2
					if sisterid in justbornsampled:
						newlin=mergeLineagesMT(lineages[i][j],justbornsampled[sisterid])
						lincur[i].append(newlin)
						del justbornsampled[sisterid]
						del justbornsampled[anc]
				else:
					extendLineage(lineages[i][j],self.dt,0)
					lincur[i].append(lineages[i][j])

			for j in justbornsampled:
				lincur[i].append(justbornsampled[j])

		return lincur
			





	def simStepForward(self,typedist):
		typedistAfterBD=self.simBDStepForward(typedist)
		typedistAfterTr=self.simTrStepForward(typedistAfterBD)
		return typedistAfterTr

	def simulateForward(self,typedist):
		self.tdStart=typedist
		tdCur=typedist
		nsteps=int((self.tEnd-self.tStart)/self.dt)
		# print nsteps

		#print 'tstart='+str(self.tStart)+'tend='+str(self.tEnd)+'nsteps='+str(nsteps)

		for i in range(nsteps):
			tdCur=self.simStepForward(tdCur)
			self.popSizeSteps.append(tdCur)
			# print tdCur

		self.tdEnd=tdCur
	
	def simStepBackward(self,step,lineages):
		#print "lins:"+str(lineages)
		linsBeforeTr=self.simTrStepBackward(step,lineages)
		#print "linsBeforeTr:"+str(linsBeforeTr)
		linsBeforeBD=self.simBDStepBackward(step,linsBeforeTr)
		#print "linsBeforeBD:"+str(linsBeforeBD)
		return linsBeforeBD

	def simulateBackward(self,lineages):

		nsteps=len(self.simBDSteps)

		# print 'BACKWARD: tstart='+str(self.tStart)+'tend='+str(self.tEnd)+'nsteps='+str(nsteps)
		
		lincur=lineages
		for i in range(nsteps-1,-1,-1):
			lincur=self.simStepBackward(i,lincur)

		return lincur


	def sampleAtBottom(self,nSamples):
		self.Nsampled=nSamples
		typedSample=sampleMultiVariateHyperGeometric(self.tdEnd,self.Nsampled)
		lineages=[ [] for i in range(self.nStates)]
		for i in range(len(typedSample)):
			for j in range(typedSample[i]):
				newnode=initLineageMT(i)
				lineages[i].append(newnode)
		self.bottomLineages=lineages
		return lineages