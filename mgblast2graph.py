#!/usr/bin/env python3

"""This is Python reimplementation of R function mgblast2graph"""

import sys
import collections


# BLASTn output format 6: http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
BlastEntry = collections.namedtuple("BlastEntry", "query subject weight "
												  "q_length q_start q_end " 
												  "s_length s_start s_end sign")




def mgblast2graph(in_blastFile, in_seqFile):
	blastEntries = loadBlastData(in_blastFile)
	




def loadBlastData(blastFile):
	result = []
	with open(blastFile, "r", encoding = "utf8") as file:
		for line in file:
			values = line[:-1].split("\t")
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





if __name__ == '__main__':
	"""Run the code on sample data"""

	params = {}
	inputFolder = "input-data/"
	params["in_blastFile"] = inputFolder + "blast.csv"
	params["in_seqFile"] = inputFolder + "reads-fas"

	mgblast2graph(**params)