#!/usr/bin/env python3

import igraph
import re

from mgblast2graph import *

# načtení základu grafu
spanningTree = igraph.Graph.Load("testFilterData/spanning.graph", format = "pajek")

# načtení jmen
with open("testFilterData/names.txt") as file:
	for num, line in enumerate(file):
		vertName = line.split("\"")[1]
		spanningTree.vs[num]["name"] = vertName

# načtení znamínek
with open("testFilterData/signs.txt") as file:
	for num, line in enumerate(file):
		sign = line.split()[1]
		spanningTree.es[num]["sign"] = int(sign)




blastEntries = loadBlastData("input-data/blast.csv")
blastEntries = filterBlastEntries(blastEntries)



reverseComplements = {int(vertex["name"]) for vertex in getNegativeEdgeVertices(spanningTree)}
# print("reverseComplements:", reverseComplements, "\n")
similarityTable, notfit = switchReversed(blastEntries, reverseComplements)

# test for similarity table correctness
df2 = []
with open("testFilterData/df2.txt") as file:
    for line in file:
        line = line.split()
        df2.append(SimilarityEntry(seq1 = int(line[1]), seq2 = int(line[2]), sign = int(line[4])))

similarityTable = set(similarityTable)
df2 = set(df2)
# print(df2 == similarityTable)



# comparison with R results
# with open("testFilterData/flipNamesResult.txt") as file:
#     fileData = file.read()

# flipNamesRegEx = re.compile(r"\"(\d+?)\"")
# flipNames = {int(x) for x in flipNamesRegEx.findall(fileData)}
# print("flipNames:", flipNames, "\n")
# print("Comparison:", reverseComplements == flipNames)


resultGraph = createResultGraph(similarityTable, notfit, reverseComplements)
clusters    = resultGraph.clusters(mode = "STRONG")
saveGraphPicture(resultGraph, createLayout(resultGraph), "result.png", "thumb.png")

# print(clusters)


