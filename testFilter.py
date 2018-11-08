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
print("reverseComplements:", reverseComplements, "\n")
similarityTable, notfit = switchReversed(blastEntries, reverseComplements)
# print("notfit:", notfit)



# comparison with R results
with open("testFilterData/flipNamesResult.txt") as file:
    fileData = file.read()

flipNamesRegEx = re.compile(r"\"(\d+?)\"")
flipNames = {int(x) for x in flipNamesRegEx.findall(fileData)}
print("flipNames:", flipNames, "\n")
print("Comparison:", reverseComplements == flipNames)