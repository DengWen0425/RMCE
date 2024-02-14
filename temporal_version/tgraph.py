import time
import copy
import bisect
import collections
from functools import cmp_to_key

class TGraph:

	def __init__(self, delta, pivot_method):
		self.begin = float('inf')
		self.end = 0
		self.lifetime = 0
		self.delta = delta
		self.o = 0
		self.pivot_method = pivot_method
		self.vertices = set()
		self.adj = {}
		self.e2intervals = {}
		self.clique_count = 0
		self.depth = 0
		self.Q1 = set()
		self.Q2 = set()

	def from_file(self, filename):
		with open(filename, "r") as f:
			lines = f.readlines()
			for line in lines:
				contents = line.split()
				t = int(contents[0])
				u = contents[1].strip()
				v = contents[2].strip()
				u, v = int(u), int(v)
				if u == v:
					continue
				if u not in self.adj:
					self.vertices.add(u)
					self.adj[u] = set()
				if v not in self.adj:
					self.vertices.add(v)
					self.adj[v] = set()
				(u, v) = (v, u) if u > v else (u, v)
				if (u, v) not in self.e2intervals:
					self.e2intervals[(u, v)] = []
				self.e2intervals[(u, v)].append(t)
				# if v not in self.adj[u]:
				self.adj[u].add(v)
				self.adj[v].add(u)
				self.begin = min(self.begin, t)
				self.end = max(self.end, t)
		
		self.begin -= self.delta
		self.end += self.delta
		self.lifetime = self.end - self.begin
		total = len(self.e2intervals.keys())
		for i,e in enumerate(list(self.e2intervals.keys())):
			u, v = e
			isD = False
			isD = self.checkCN(u,v)
			isU = len(self.adj[u]) == 1
			isV = len(self.adj[v]) == 1
			# isU, isV = False, False
			intervals = []
			eList = sorted(self.e2intervals[e])
			start = max(self.begin, eList.pop(0) - self.delta)
			end = min(self.begin + self.lifetime, start + self.delta * 2)
			for t in eList:
				if (end+1) < t:
					if isU or isV or isD:
						self.o += 1
					intervals.append((start,end))
					start = max(self.begin, t - self.delta)
					end = min(self.begin + self.lifetime,t + self.delta)
				else:
					end = min(self.begin + self.lifetime, t + self.delta)
			intervals.append((start,end))
			if not isU and not isV and not isD:
				self.e2intervals[e] = intervals
			else:
				self.o += 1
				self.adj[v].remove(u)
				self.adj[u].remove(v)
				if isU:
					self.vertices.remove(u)
				if isV:
					self.vertices.remove(v)
				del self.e2intervals[e]
		for v in self.vertices:
			if len(self.adj[v]) == 1:
				self.Q1.add(v)
			if len(self.adj[v]) == 2:
				self.Q2.add(v)
		return 

	# initial call
	def run(self):
		self.globalReduction()
		R = ([], (self.begin, self.begin + self.lifetime))
		P = dict()
		X = dict()
		for v in self.vertices:
			P[v] = [(self.begin, self.begin + self.lifetime)]
		self.recur(P, R, X)


	def recur(self, P, R, X):
		self.depth += 1
		if not P and not X and self.depth > 2:
			self.clique_count += 1
             
		dismiss = []
		if self.pivot_method > 0:
			recompute = True
			while recompute:
				pivots, pivotIntervals, recompute = self.computePivotInformationWithReduction(R,P,X, copy.deepcopy(P))
			if len(pivots) > 1:
				if self.pivot_method == 1:
					dismiss = pivotingOne(pivots, pivotIntervals)
				elif self.pivot_method == 2:
					dismiss = pivotingOneMax(pivots, pivotIntervals)
				elif self.pivot_method == 3:
					dismiss = pivotingMany(pivots, pivotIntervals)
				elif self.pivot_method == 4:
					dismiss = pivotingGreedy(pivots, pivotIntervals)
				elif self.pivot_method == 5:
					dismiss = pivotingOptimal(pivots, pivotIntervals)
			elif len(pivots) == 1:
				_, dismiss = pivots.popitem()

		P_support = copy.deepcopy(P)
		while len(P_support) != 0:
			v, In = P_support.popitem()

			if v not in X: X[v] = []

			while len(In) != 0:
				interval = In.pop()
				if (v, interval) not in dismiss:
					R_new = None
					P_new = self.cut(v, *interval, copy.deepcopy(P))
					X_new = self.cut(v, *interval, copy.deepcopy(X))
					self.recur(P_new, R_new, X_new)
					self.depth -= 1
					P[v].remove(interval)
					bisect.insort(X[v], interval)
			if len(P[v]) == 0: del P[v]
			if len(X[v]) == 0: del X[v]

	def computePivotInformationWithReduction(self, R, oriP, X, P):
		pivot = {}
		pivotIntervals = []
		pCnt = sum([len(P[l]) for l in P])
		recompute = False
		for v in P:
			for I in P[v]:
				can0prune = True
				pivot[(v,I)] = set()
				for w in P:
					if w != v:
						tmp, pruneCheck = self.checkIn(w, P[w], v, *I)
						if not pruneCheck:
							can0prune = False
						pivot[(v,I)].update(tmp)
				if len(pivot[(v,I)]) == 0: 
					del pivot[(v,I)]
					if can0prune:
						if self.checkX(X,v, I):
							self.clique_count += 1
						oriP[v].remove(I)
				elif len(pivot[(v,I)]) == pCnt-1: # deg-p-1
					recompute = True
					oriP[v].remove(I)
					if len(oriP[v]) == 0:
						del oriP[v]
					X = self.cut(v, *I, X)
					return pivot, pivotIntervals, recompute
				else:
					pivotIntervals.append(DeltaInterval(v,*I, pivot[(v,I)]))  
		return pivot, pivotIntervals, recompute

	def checkIn(self, w, list1, v, start, end):
		pivot = set()
		key = (w, v) if v > w else (v, w)
		if key in self.e2intervals:
			pruneCheck = True
			list2 = self.e2intervals[key]
			list1_curr = 0
			list2_curr = 0
			while list1_curr < len(list1) and list2_curr < len(list2):
				(i1_start, i1_end) = list1[list1_curr]
				if end < i1_end or start > i1_start:
					pruneCheck = False
					break
				elif end < (i1_start + self.delta):
					break
				elif i1_end < (start + self.delta):
					list1_curr += 1
					continue
				else:
					while list2_curr < len(list2):
						(i2_start, i2_end) = list2[list2_curr]
						if i2_end < (start + self.delta):
							list2_curr += 1
							continue
						if end < (i2_start + self.delta):
							list2_curr = len(list2)
							break
						if i2_end < (i1_start + self.delta):
							list2_curr += 1
							continue
						if i2_start >= i1_end: 
							list1_curr += 1
							break	
						if i1_start >= i2_start and i1_end <= i2_end:
							pivot.add((w,(i1_start, i1_end)))
							list1_curr += 1
							pruneCheck = False
							break
						if pruneCheck:
							if i2_start < i1_start and i2_end > i1_start and i2_end <= i1_end:
								a = max(i1_start, start)
								b = min(i2_end, end)
								if b - a >= self.delta:
									pruneCheck = False

							elif i2_start < i1_start and i2_end > i1_end: # dominated
								a = max(i1_start, start)
								b = min(i1_end, end)
								if b - a >= self.delta:
									pruneCheck = False

							elif i2_start >= i1_start and i2_end <= i1_end:
								a = max(i2_start, start)
								b = min(i2_end, end)
								if b - a >= self.delta:
									pruneCheck = False					

							elif i2_start >= i1_start and i2_start < i1_end and i2_end > i1_end:
								a = max(i2_start, start)
								b = min(i1_end, end)
								if b - a >= self.delta:
									pruneCheck = False

						if i1_start < i2_start and i1_end > i2_end:
							list1_curr += 1
							list2_curr += 1
							break
						if i1_start < i2_start:
							list1_curr += 1
							break
						if i1_end > i2_end:
							list2_curr += 1
							continue
		else:
			pruneCheck = True
		return pivot, pruneCheck 

	def cut(self, v, start, end, P):
		to_iterate = list(P.keys())
		for w in to_iterate:
			key = (w, v) if v > w else (v, w)
			if key in self.e2intervals:
				P[w] = self.deltaIntersect(start, end, P[w], self.e2intervals[key])
				if len(P[w]) == 0:
					del P[w]		
			else: 
				del P[w]
		return P  

	def deltaIntersect(self, start, end, list1, list2):
		list_new = []
		list1_curr = 0
		list2_curr = 0

		while list1_curr < len(list1) and list2_curr < len(list2):
			(i1_start, i1_end) = list1[list1_curr]
			if end < (i1_start + self.delta):
				break
			if i1_end < (start + self.delta):
				list1_curr += 1
				continue
			else:
				while list2_curr < len(list2):
					(i2_start, i2_end) = list2[list2_curr]						
					if i2_end < (start + self.delta):
						list2_curr += 1
						continue
					if end < (i2_start + self.delta):
						list2_curr = len(list2)
						break
					if i2_end < (i1_start + self.delta):
						list2_curr += 1
						continue
					if i2_start < i1_start and i2_end > i1_start and i2_end <= i1_end:
						a = max(i1_start, start)
						b = min(i2_end, end)
						if b - a >= self.delta:
							list_new.append((a,b)) 
						list2_curr += 1
						continue

					if i2_start < i1_start and i2_end > i1_end: # dominated
						a = max(i1_start, start)
						b = min(i1_end, end)
						if b - a >= self.delta:
							list_new.append((a,b)) 
						list1_curr += 1
						break

					if i2_start >= i1_start and i2_end <= i1_end:
						a = max(i2_start, start)
						b = min(i2_end, end)
						if b - a >= self.delta:
							list_new.append((a,b))
						list2_curr += 1
						continue						

					if i2_start >= i1_start and i2_start < i1_end and i2_end > i1_end:
						a = max(i2_start, start)
						b = min(i1_end, end)
						if b - a >= self.delta:
							list_new.append((a,b))
						list1_curr += 1
						break

					if i2_start >= i1_end: 
						list1_curr += 1
						break		

		return list_new 

	def checkCN(self, u, v):
		for w in self.adj[u]:
			if w in self.adj[v]:
				return False
		return True

	def checkX(self, X, v, I):
		for u in X.keys():
			if v in self.adj[u]:
				key = (u, v) if v > u else (v, u)
				if key in self.e2intervals:
					for I2 in self.e2intervals[key]:
						if I[0] >= I2[0] and I[1] <= I2[1]:
							return False
		return True


	def globalReduction(self):
		while len(self.Q1) != 0 or len(self.Q2) != 0:
			while len(self.Q1) != 0:
				v = self.Q1.pop()
				if v not in self.vertices:
					continue
				# u = self.adj[v][0]
				u = self.adj[v].pop()
				key = (u, v) if v > u else (v, u)
				self.clique_count += len(self.e2intervals[key])
				del self.e2intervals[key]
				self.vertices.remove(v)
				self.adj[u].remove(v)
				self.adj.pop(v)
				if len(self.adj[u]) == 1:
					self.Q2.remove(u)
					self.Q1.add(u)
				elif len(self.adj[u]) == 0:
					self.vertices.remove(u)
					self.adj.pop(u)
				elif len(self.adj[u]) == 2:
					self.Q2.add(u)
			while len(self.Q2) != 0:
				v = self.Q2.pop()
				u,w = list(self.adj[v])
				if u not in self.adj[w]:
					key1 = (u, v) if v > u else (v, u)
					key2 = (w, v) if v > w else (v, w)
					self.clique_count += len(self.e2intervals[key1])
					self.clique_count += len(self.e2intervals[key2])
					del self.e2intervals[key1]
					del self.e2intervals[key2]
					self.vertices.remove(v)
					self.adj[u].remove(v)
					self.adj[w].remove(v)
					self.adj.pop(v)
					if len(self.adj[u]) == 1:
						self.Q2.remove(u)
						self.Q1.add(u)
					elif len(self.adj[u]) == 0:
						self.vertices.remove(u)
						self.adj.pop(u)
					elif len(self.adj[u]) == 2:
						self.Q2.add(u)
					if len(self.adj[w]) == 1:
						self.Q2.remove(w)
						self.Q1.add(w)
					elif len(self.adj[w]) == 0:
						self.vertices.remove(w)
						self.adj.pop(w)
					elif len(self.adj[w]) == 2:
						self.Q2.add(w)
		return


class DeltaInterval:
	def __init__(self, v,s,e,d):
		self.vertex = v
		self.startTime = s
		self.finishTime = e
		self.dismissElements = d

def pivotingOne(pivots, pivotIntervals):
	_, dismissElements = pivots.popitem()
	return dismissElements

def pivotingOneMax(pivots, pivotIntervals):
	maxValue = 0
	maxElement = 0
	dismissElements = set()
	for el in pivots:
		l = len(pivots[el])
		if l  > maxValue:
			maxValue = l
			maxElement = el
	if maxValue > 0:
		dismissElements = pivots[maxElement]
	return dismissElements

def pivotingMany(pivots, pivotIntervals):

	pivot_elements = []
	dismissElements = set()
	n_chosen = 0
	while len(pivots) > 0:
		act_el, new_dismiss_elements = pivots.popitem()

		_, (ae_s, ae_e) = act_el
		
		eligible = True

		for piv_el in pivot_elements:
			_, (pe_s, pe_e) = piv_el
			a = max(ae_s, pe_s)
			b = min(ae_e, pe_e)
			if a <= b:
					eligible = False
					break

		if eligible:
			pivot_elements.append(act_el)
			dismissElements = dismissElements.union(new_dismiss_elements)
			n_chosen += 1
	return dismissElements

def find_max_el(pivots, pivotIntervals):
    
	maxValue = 0
	maxElement = 0
	dismissElements = set()
	for el in pivots:
		l = len(pivots[el])
		if l  > maxValue:
			maxValue = l
			maxElement = el

	dismissElements = pivots[maxElement]

	return maxElement, dismissElements

def pivotingGreedy(pivots, pivotIntervals):
	pivot_elements = []
	dismiss_el = set()
	n_chosen = 0
	
	while len(pivots) > 0:
		max_el, new_dismiss_el = find_max_el(pivots, pivotIntervals)
		n_chosen += 1

		dismiss_el = dismiss_el.union(new_dismiss_el)

		_, (me_s, me_e) = max_el

		rem_pivots = []
		for piv_el in pivots:
			_, (pe_s, pe_e) = piv_el
			a = max(me_s, pe_s)
			b = min(me_e, pe_e)
			if a <= b:
				rem_pivots = rem_pivots + [piv_el]

		for re in rem_pivots:
			pivots.pop(re)
	return dismiss_el

def pivotingOptimal(pivots, pivotIntervals):
	n_ivals = len(pivotIntervals)
	def interval_compare(x, y):
		if x.finishTime == y.finishTime:
			if x.vertex <= y.vertex:
				return -1
			else:
				return 1
		else:
			return x.finishTime - y.finishTime

	pivotIntervals.sort(key = cmp_to_key(interval_compare))

	succinct_start = []

	for i in range(0, n_ivals):
		item = (pivotIntervals[i], i)
		succinct_start.append(item)
			
	succinct_start.sort(key = lambda item: item[0].startTime)

	previous = [-1] * n_ivals
	curr = 0
	for i in range(1, n_ivals):
		while pivotIntervals[curr].finishTime < succinct_start[i][0].startTime:
			curr += 1
		if pivotIntervals[curr].startTime == pivotIntervals[curr].finishTime:
			previous[succinct_start[i][1]] = curr - 1
		elif pivotIntervals[curr].finishTime > succinct_start[i][0].startTime:
			previous[succinct_start[i][1]] = curr - 1
		else:
			previous[succinct_start[i][1]] = curr
	opt = [0] * n_ivals
	opt[0] = len(pivotIntervals[0].dismissElements)
	chosen = [-1] * n_ivals
	chosen[0] = (True, -1)

	for i in range(1, n_ivals):
		n_dismiss = len(pivotIntervals[i].dismissElements)
		if previous[i] < 0:
			if opt[i - 1] <= n_dismiss:
				opt[i] = n_dismiss
				chosen[i] = (True, -1)
			else:
				opt[i] = opt[i - 1]
				chosen[i] = (False, i - 1)
		else:
			if opt[i - 1] <= opt[previous[i]] + n_dismiss:
				opt[i] = opt[previous[i]] + n_dismiss
				chosen[i] = (True, previous[i])
			else:
				opt[i] = opt[i - 1]
				chosen[i] = (False, i - 1)
	dismiss_elems = set()
	i = n_ivals - 1
	n_chosen = 0
	
	while i >= 0:
		ival_chosen, _ = chosen[i]
		if ival_chosen:
			ival = pivotIntervals[i]
			dismiss_elems.update(ival.dismissElements)
			i = previous[i]
			n_chosen += 1
		else:
			i -= 1
	return dismiss_elems