#!/usr/bin/env python3

"""This is Python reimplementation of R function mgblast2graph"""

import sys
import collections

## TYPE DEFINITIONS ##

# BLASTn output format 6: http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
BlastEntry = collections.namedtuple("BlastEntry", "query subject weight "
												  "q_length q_start q_end " 
												  "s_length s_start s_end sign")

Sequence = collections.namedtuple("Sequence", "description sequence")

class FileError(Exception): pass


## FUNCTIONS ##

def mgblast2graph(blastFileName, seqFileName):
	blastEntries = loadBlastData(blastFileName)
	blastEntries = filterBlastEntries(blastEntries)
	sequences = loadSequences(seqFileName)



def loadBlastData(blastFileName):
	"""Loads data from blast file"""

	result = []
	with open(blastFileName, "r", encoding = "utf8") as file:
		for line in file:
			if line[-1] == "\n":
				line = line[:-1]
			values = line.split("\t")
			entry = BlastEntry(
					query    = int(values[0]),
					subject  = int(values[1]),
					weight   = float(values[2]),

					q_length = int(values[3]),
					q_start  = int(values[4]),
					q_end    = int(values[5]),

					s_length = int(values[6]),
					s_start  = int(values[7]),
					s_end    = int(values[8]),

					sign     = values[11]
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
			if letters - {"A", "C", "G", "T"}:
				raise FileError("Incorrect sequence data")

			result.append(Sequence(description = line1[1:], sequence = line2))

	return result





if __name__ == '__main__':
	"""Run the code on sample data"""

	#doctesting
	import doctest
	res = doctest.testmod()
	if res.failed != 0:
		sys.exit(1)

	#set params
	params = {}
	inputFolder = "input-data/"
	params["blastFileName"] = inputFolder + "blast.csv"
	params["seqFileName"] = inputFolder + "reads.fas"

	mgblast2graph(**params)
	