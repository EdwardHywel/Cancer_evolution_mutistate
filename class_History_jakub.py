#!/usr/bin/env python

class History:

	tree=None # each node corresponds to a deme

	def __init__(self,tree):
		self.tree=dendropy.Tree(tree)

	def dividePopBetweenChildren(self,node):

		if node.is_leaf():
			return #no need to do anything

		children=[child for child in node.child_nodes()]
		totalweight=sum(map(lambda x:x.deme.weight,children))
		relweights=map(lambda x: x.deme.weight/totalweight,children)

		if isinstance(node.deme,MultiTypeDeme):# for multitype, we have to keep track of how many from each type are included in the descendant populations
			for state in range(node.deme.nStates):
				division=list(np.random.multinomial(node.deme.tdEnd[state],relweights))
				for i in range(len(children)):
					children[i].deme.tdStart[state]=division[i]
		else:	
			division=list(np.random.multinomial(node.deme.Nend,relweights))
			for i in range(len(children)):
				children[i].deme.Nstart=division[i]


	def simulateForwardOld(self,Nfounders):

		for node in self.tree.preorder_node_iter():
			if node.parent_node==None:
				node.deme.Nstart=Nfounders
				node.deme.simulateForward(Nfounders)
				self.dividePopBetweenChildren(node)
				continue
			node.deme.simulateForward(node.deme.Nstart)
			self.dividePopBetweenChildren(node)

	def simulateForwardMultiType(self,Nfounders):#Nfounders is a typedist, not a num

		for node in self.tree.preorder_node_iter():
			if node.parent_node==None:
				node.deme.tdStart=Nfounders
				node.deme.simulateForward(Nfounders)
				self.dividePopBetweenChildren(node)
				continue
			node.deme.simulateForward(node.deme.tdStart)
			self.dividePopBetweenChildren(node)

	def simulateForward(self,Nfounders):
		if isinstance(self.tree.seed_node.deme,MultiTypeDeme):
			self.simulateForwardMultiType(Nfounders)
		else:
			self.simulateForwardOld(Nfounders)

	#assumes deme.Nsampled stores the required number of samples in each leaf deme
	def sampleAtBottom(self):
		for leaf in self.tree.leaf_nodes():
			leaf.deme.sampleAtBottom(leaf.deme.Nsampled)




	def simulateBackward(self):

		for node in self.tree.postorder_node_iter():
			lins=node.deme.simulateBackward(node.deme.bottomLineages)
			#print "LINEAGES:\n"
			#printLineages(lins)
			if node.parent_node!=None:
				node.parent_node.deme.bottomLineages.extend(lins)
			#print node.deme.Nend

		return lins


	def simulate(self,Nfounders):
		self.simulateForward(Nfounders)
		self.sampleAtBottom()
		lins=self.simulateBackward()
		return lins