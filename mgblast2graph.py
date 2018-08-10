#!/usr/bin/env python3

"""This is Python reimplementation of R function mgblast2graph"""

def mgblast2graph(in_blastFile, in_seqFile):
	pass


if __name__ == '__main__':
	"""Run the code on sample data"""

	params = {}
	inputFolder = "input-data/"
	params["in_blastFile"] = inputFolder + "in-blast.csv"
	params["in_seqFile"] = inputFolder + "in-reads-fas"

	mgblast2graph(**params)