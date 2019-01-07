#!/usr/bin/env python3

"""This is Python reimplementation of R function mgblast2graph"""

import sys
import collections
import math
import random
import copy
import os
import os.path

import igraph
from PIL import Image	# Pillow library


## TYPE DEFINITIONS ##

# BLASTn output format 6: http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
BlastEntry = collections.namedtuple("BlastEntry", "query subject weight "
												  "q_length q_start q_end " 
												  "s_length s_start s_end sign")

SimilarityEntry = collections.namedtuple("SimilarityEntry", "seq1 seq2 sign")

Sequence = collections.namedtuple("Sequence", "description sequence")

class FileError(Exception): pass


## FUNCTIONS ##

def mgblast2graph(blastFileName: str, seqFileName: str,
				  maxSampleVertices: int, maxSampleEdges: int,
				  paired: bool, satelliteModelFile: str,
				  pictureName: str, thumbnailName: str, outputSeqFileName: str):

	createFoldersForFiles(pictureName, thumbnailName, outputSeqFileName)

	blastEntries = loadBlastData(blastFileName)
	blastEntries = filterBlastEntries(blastEntries)
	sequences    = loadSequences(seqFileName)

	# compute pair complettness index before sampling
	pairCompletenessIndex    = getPairCompletenessIndex(sequences) if paired else 0
	sequences, blastEntries = createSample(sequences, blastEntries,
										   maxSampleVertices, maxSampleEdges)

	graph  = createGraph(sequences, blastEntries)
	layout = createLayout(graph)
	# original script saves graph and layout to *.GL file here
	# ommiting for now
	saveGraphPicture(graph, layout, pictureName, thumbnailName)

	if not graph.is_connected():	# sequences are not filtered?
		graph, blastEntries = getLargestComponent(graph, blastEntries)

	spanningTree = graph.spanning_tree(weights = [edge["weight"] for edge in graph.es])
	# original script tries to make alternative spanning trees here
	# in case that "suboptimal solution is found", ignoring for now

	reverseComplements = {vertex["name"] for vertex in getNegativeEdgeVertices(spanningTree)}
	similarityTable, notfit = switchReversed(blastEntries, reverseComplements)

	resultGraph = createResultGraph(similarityTable, notfit, reverseComplements)
	clusters    = resultGraph.clusters(mode = "STRONG")
	membership  = vertexToClusterMembership(clusters)

	resultSequences = alterSequences(sequences, reverseComplements, notfit)
	saveSequencesAndClusterData(resultSequences, resultGraph, membership, outputSeqFileName)

	# GRAPH INFO COMPUTATION
	# escore is sum of entries with sign 1 divided by all entries
	escore = sum(entry.sign for entry in similarityTable if entry.sign == 1)/len(similarityTable)
	coverage = len(resultSequences)/len(sequences)
	loopIndex = max([len(cluster) for cluster in clusters])/len(resultGraph.vs)

	# calculate satellite probability
	matrix, cutoff = loadSatelliteModel(satelliteModelFile)
	sattProb = getSattProbability(matrix, loopIndex, pairCompletenessIndex)
		# warning: isSatt has values "Putative Satellite"/"" not True/False

	graphInfo = {
		"escore"                : escore,
		"escore_mts"            : escoreMstComputation(spanningTree),
		"coverage"              : coverage,
		"loop_index"            : loopIndex,
		"pair_completeness"     : pairCompletenessIndex,
		"graph_file"            : None,	# ignoring for now
		"oriented_sequences"    : outputSeqFileName,
		"vcount"                : len(resultGraph.vs),
		"ecount"                : len(resultGraph.es),
		"satellite_probability" : sattProb,
		"satellite"             : "Putative Satellite" if sattProb > cutoff else ""
	}

	return graphInfo



def createFoldersForFiles(*filePaths):
	for path in set(map(os.path.dirname, filePaths)):
		if os.path.exists(path): continue
		os.makedirs(path)
		

def loadBlastData(blastFileName):
	"""Loads data from blast file"""

	result = []
	with open(blastFileName, "r", encoding = "utf8") as file:
		for line in file:
			if line[-1] == "\n":
				line = line[:-1]
			values = line.split("\t")

			if values[11] not in {"+", "-"}:
				raise FileError("Incorrect sign in blast file")

			entry = BlastEntry(
					query    = values[0],		# ID of sequence has to be string in order to 
					subject  = values[1],		# allow for paired information
					weight   = float(values[2]),

					q_length = int(values[3]),
					q_start  = int(values[4]),
					q_end    = int(values[5]),

					s_length = int(values[6]),
					s_start  = int(values[7]),
					s_end    = int(values[8]),

					sign     = +1 if values[11] == "+" else -1
				)
			result.append(entry)
	return result


def filterBlastEntries(blastEntries):
	"""Filteres duplicated entries.
	
	Since both query and subject references same set of sequences,
	same number is considered same sequence no matter whether it is saved in
	query or subject field.

	>>> entries = [
	... BlastEntry(query = 1, subject = 2, weight = None, q_length = None, 
	...			   q_start = None, q_end = None, s_length = None, s_start = None, 
	...			   s_end = None, sign = None),
	... BlastEntry(query = 2, subject = 1, weight = None, q_length = None, 
	...			   q_start = None, q_end = None, s_length = None, s_start = None, 
	...			   s_end = None, sign = None),
	... BlastEntry(query = 1, subject = 2, weight = None, q_length = None, 
	...			   q_start = None, q_end = None, s_length = None, s_start = None, 
	...			   s_end = None, sign = None)]

	>>> entries = filterBlastEntries(entries)
	>>> len(entries)
	1
	>>> entries[0][:2]
	(1, 2)

	"""

	uniqueKeys = set()
	filteredEntries = []

	for entry in blastEntries:
		keys = {(entry.query, entry.subject), (entry.subject, entry.query)}
		if not (keys & uniqueKeys):
			uniqueKeys.add((entry.query, entry.subject))
			filteredEntries.append(entry)

	return filteredEntries


def loadSequences(seqFileName):
	"""Loads sequences from FASTA file"""

	result = []
	with open(seqFileName, "r", encoding = "utf8") as file:
		while True:
			#get annotation
			line1 = file.readline().strip()
			if len(line1) == 0: break
			if line1[0] != ">" or len(line1) == 1:
				raise FileError("Broken file structure")

			#get sequence data
			line2 = file.readline().strip()
			if len(line2) == 0:
				raise FileError("Broken file structure")

			#ACGT test
			letters = set(line2)
			if letters - {"A", "C", "G", "T", "N"}:
				raise FileError("Incorrect sequence data")

			result.append(Sequence(description = line1[1:], sequence = line2))

	return result


def getPairCompletenessIndex(sequences):
	"""Compute portion of broken read pairs

	P = Nc/(Nc + Ni)
	P - pair completness index
	Nc - number of complete read pairs
	Ni - number of broken pairs
	"""

	# first we count occurence of the same sequence description
	# while ignoring last digit of the description, I don't know why it was
	# like that in the original script
	nameCount = collections.defaultdict(lambda: 0)
	for sequence in sequences:
		nameCount[sequence.description[:-1]] += 1

	# second, we count occurences of the amounts computed in previous step
	# like: how many times we get 1 occurence, how many 2 occurences ...
	occurenceCount = collections.defaultdict(lambda: 0)
	for counts in nameCount.values():
		occurenceCount[counts] += 1

	# the result is 1 - number of occurences that has happenec only once
	# divided by sum of all occurences, that has happened any amount of times
	return 1 - occurenceCount[1]/sum(occurenceCount.values())


def createSample(sequences, blastEntries, maxVertices, maxEdges):
	"""Samples input data if they are too large"""

	vertices = len(sequences)
	edges = len(blastEntries)

	#check whether sampling is necessary
	if maxVertices > vertices and maxEdges > edges:
		return sequences, blastEntries

	sampleSize = estimateSampleSize(vertices, edges, maxVertices, maxEdges)
	sequences = random.sample(sequences, sampleSize)

	filteredEntries = []
	seqIDs = {seq.description for seq in sequences}
	for blastEntry in blastEntries:
		if blastEntry.query in seqIDs and blastEntry.subject in seqIDs:
			filteredEntries.append(blastEntry)

	return sequences, filteredEntries


def estimateSampleSize(vertices, edges, maxVertices, maxEdges):
	"""Estimates sample size for graph"""

	# I think this whole function ca be simplified just to:
	#	result = vertices if vertices <= maxVertices else maxVertices
	d = (2 * edges)/(vertices * (vertices - 1))
	eEst = (maxVertices * (maxVertices - 1) * d)/2
	nEst = (d + math.sqrt((d ** 2) + 8 * d * maxEdges))/(2 * d)

	if eEst >= maxEdges:
		N = round(nEst)
		E = round((N * (N - 1) * d)/2)

	if nEst >= maxVertices:
		N = maxVertices
		E = round((N * (N - 1) * d)/2)

	return N


def createGraph(sequences, blastEntries):
	"""Creates graph where vertices are sequences (their names)
	and edges are created from blast entries"""

	graph = igraph.Graph(directed = False)

	for seq in sequences:
		graph.add_vertex(name = seq.description)

	for entry in blastEntries:
		graph.add_edge(entry.query, entry.subject, 
					   weight = entry.weight, sign = entry.sign)

	return graph


def createLayout(graph):
	if graph.ecount() < 2000000:
		# missing OGDF computation
		# original script methods.R line 921, 922
		# ignored for now, because of the need to create interface for C++
		# ogdf library
		layout = graph.layout_fruchterman_reingold()
	else:
		layout = graph.layout_fruchterman_reingold()

	return layout


def saveGraphPicture(graph, layout, pictureName, thumbnailName):
	"""Saves graph as a picture with it's thumbnail"""

	# Plotting has to be setup:
	# 1) install cairo library (using standard package manager)
	# 2) python interface for cairo - cairocffi: pip install cairocffi

	if ".png" != pictureName[-4:]:
		raise ValueError("Picture must be PNG")
	if ".png" != thumbnailName[-4:]:
		raise ValueError("Thumbnail must be PNG")

	igraph.plot(graph,
				layout = layout,
				target = pictureName,
				bbox = (900, 900), 
				vertex_size = 2,
				vertex_shape = "circle",
				vertex_color = "#000000A0",
				edge_color = "#00000040",
				edge_curved = False,
				edge_width = 1,
				vertex_label = None,
				vertex_frame_color = None)

	with Image.open(pictureName) as image:
		thumbnail = image.resize((100, 100), Image.ANTIALIAS)
		thumbnail.save(thumbnailName)


def getLargestComponent(graph, blastEntries):
	clusters = graph.clusters()
	biggestSubgraph = clusters.giant()
	
	filteredEntries = []
	vertexNames = {int(vertex["name"]) for vertex in biggestSubgraph.vs}

	for entry in blastEntries:
		if entry.query in vertexNames and entry.subject in vertexNames:
			filteredEntries.append(entry)

	return biggestSubgraph, filteredEntries
	

def getNegativeEdgeVertices(spanningTree):
	"""Returns a list of vertices, that are connected to a predecessor vertex
	with a negative edge.
	"""
	
	result = []
	spanningTree = spanningTree.copy()

	for vertex, parent in depthFirstSearch(spanningTree, 0):
		if vertex == parent: continue	# ignoring root of the search

		edge = spanningTree.es.select(_between = ((vertex,),(parent,)))[0]
		if edge["sign"] == -1:
			result.append(spanningTree.vs[vertex])

			# because the resulting vertex would be switched, switching all 
			# edges connected to it
			edges = list(spanningTree.es.select(_source = vertex))
			edges += list(spanningTree.es.select(_target = vertex))
			for edge in edges: edge["sign"] *= -1
		
	return result


def escoreMstComputation(spanningTree):
	"""Returns sum of all altered edges divided by their amount.
	Not sure what mst means in R version.

	This function uses the same algotithm as getNegativeEdgeVertices
	to create altered spanning tree.
	This of course is bad design to repeat similar algorithm twice and
	should be redone in next iteration.
	"""
	
	spanningTree = spanningTree.copy()

	for vertex, parent in depthFirstSearch(spanningTree, 0):
		if vertex == parent: continue	# ignoring root of the search

		edge = spanningTree.es.select(_between = ((vertex,),(parent,)))[0]
		if edge["sign"] == -1:

			# because the resulting vertex would be switched, switching all 
			# edges connected to it
			edges = list(spanningTree.es.select(_source = vertex))
			edges += list(spanningTree.es.select(_target = vertex))
			for edge in edges: edge["sign"] *= -1
		
	return sum([edge["sign"] for edge in spanningTree.es])/len(spanningTree.es)


def depthFirstSearch(graph, startVertexNumber):
	"""The depth first search algorithm.
	Returns list of tuples - (vertexNumber, parentVertexNumber)
	"""

	stack = [(startVertexNumber, startVertexNumber)]
	visitedVertices = set()

	while stack:
		vertexNum, parentNum = stack.pop(0)
		if vertexNum not in visitedVertices:
			visitedVertices.add(vertexNum)
			yield vertexNum, parentNum

			vertex = graph.vs[vertexNum]
			neighbors = [(neighbor.index, vertexNum) for neighbor in vertex.neighbors()]
			stack = neighbors + stack


def switchReversed(blastEntries, reverseComplements):
	similarityTable = []
	notfit = set()

	for entry in blastEntries:
		sign = entry.sign

		queryMean = (entry.q_start + entry.q_end)/2
		if entry.query in reverseComplements:
			queryMean = entry.q_length - queryMean
			sign *= -1

		subjectMean = (entry.s_start + entry.s_end)/2
		if entry.subject in reverseComplements:
			subjectMean = entry.s_length - subjectMean
			sign *= -1

		# set new order of sequences in entry
		if subjectMean > queryMean: 
			seq1, seq2 = entry.subject, entry.query
		else: 
			seq1, seq2 = entry.query, entry.subject
		
		similarityTable.append(SimilarityEntry(seq1, seq2, sign))
		if sign == -1: notfit |= {seq1, seq2}

	return similarityTable, notfit


def createResultGraph(similarityTable, notfit, reverseComplements):
	"""Creates graph from vertices, that are not in nofit
	and edges among those vertices
	"""

	graph = igraph.Graph(directed = True)
	vertices = {entry.seq1 for entry in similarityTable} | \
			   {entry.seq2 for entry in similarityTable}
	vertices -= notfit
	
	for vertex in vertices:
		graph.add_vertex(name = vertex, complement = vertex in reverseComplements)

	for entry in similarityTable:
		if entry.seq1 not in notfit and entry.seq2 not in notfit:
			graph.add_edge(entry.seq1, entry.seq2, sign = entry.sign)
	
	return graph


def alterSequences(sequences, reverseComplements, notfit):
	"""Removes nofit sequences and reverses reverse complements"""

	result = []
	complMap = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

	for seq in sequences:
		seqID = seq.description
		seqData = seq.sequence

		if seqID in notfit: continue
		if seqID in reverseComplements:
			bases = list(seqData)
			bases = reversed([complMap[base] for base in bases])
			seqData = "".join(bases)

		result.append(Sequence(seqID, seqData))
	return result


def vertexToClusterMembership(clusters):
	"""Returns dictionary with mapping: vertexIndex -> clusterIndex"""

	membership = {}

	for index, vertices in enumerate(clusters, 1):
		for vertex in vertices:
			membership[vertex] = index

	return membership


def saveSequencesAndClusterData(sequences, resultGraph, membership, fileName):
	"""Saves sequences with info about cluster they belong to"""

	# creating sort map - sorts clusters by size
	# this is solely for mimicking the exact result of original R version
	# that saves biggest clusters with smallest cluster indexes
	sortMap = {}
	counter = collections.defaultdict(int)	# counting how many vertices each cluster has
	for vertex in membership:
		counter[membership[vertex]] += 1

	clusterSizes = sorted(counter.items(), key = lambda x: x[1], reverse = True)
	for newIndex, clusterInfo in enumerate(clusterSizes, start = 1):
		sortMap[clusterInfo[0]] = newIndex


	with open(fileName, "w", encoding = "utf8") as file:
		for seq in sequences:	# does sequences have to be filtered?
			vertexIndex = resultGraph.vs.find(seq.description).index
			clusterIndex = sortMap[membership[vertexIndex]]
			file.write("".join((">", seq.description, " ", str(clusterIndex), "\n")))
			file.write(seq.sequence + "\n")


def loadSatelliteModel(satelliteModelFile):
	with open(satelliteModelFile, "r", encoding = "utf8") as file:
		# get cutoff
		line = file.readline()
		if "cutoff" not in line[:6]:
			raise ValueError("Incorrect data in satelliteModelFile")
		cutoff = float(line.split("=")[1])
		
		# ignore description in file
		for line in file:
			if "probability_matrix" in line: break

		# get probability matrix
		matrix = {}
		for lineIndex, line in enumerate(file, start = 1):
			line = line.split()
			for columnIndex, value in enumerate(line, start = 1):
				matrix[lineIndex, columnIndex] = float(value)

	return matrix, cutoff


def getSattProbability(matrix, loopIndex, pairCompletenessIndex):
	x, y = loopIndex, pairCompletenessIndex

	# sattelite probability
	N = int(math.sqrt(len(matrix)))	# number of lines or columns
	i = round(x * (N - 1)) + 1		# rescaling params for matrix
	j = round(y * (N - 1)) + 1
	probability = matrix[i, j]

	return probability



if __name__ == '__main__':
	"""Run the code for testing purposes on sample data"""

	#doctesting
	import doctest
	res = doctest.testmod()
	if res.failed != 0:
		sys.exit(1)

	#set params
	params = {}
	inputFolder  = "input-data/"
	outputFolder = "output-data/"

	params["blastFileName"]      = inputFolder + "blastPaired.csv"
	params["seqFileName"]        = inputFolder + "readsPaired.fas"
	params["satelliteModelFile"] = inputFolder + "satellite_model.txt"
	params["maxSampleVertices"]  = 40000
	params["maxSampleEdges"]     = 20000000
	params["paired"]             = True
	params["pictureName"]        = outputFolder + "graphPicture.png"
	params["thumbnailName"]      = outputFolder + "graphThumbnail.png"
	params["outputSeqFileName"]  = outputFolder + "sequences.fasta"

	result = mgblast2graph(**params)

	print()
	for key in sorted(result.keys()):
		print(key, ":", result[key])
	print()
