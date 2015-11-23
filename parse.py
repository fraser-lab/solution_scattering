""" This part of the package is for loading data of various types and then 
making traces

Benjamin Barad
"""
from numpy import recarray, load
from trace import Trace


def parse_tpkl(filename):
	"""Loads tpkl files and generates a corresponding Trace object. Requires
	table.py from the Anfinrud lab, which we will not distribute.
	"""
	from table import table
	data = load(filename)
	q = data['q']
	S = data['S']
	SA = data['SA']
	return Trace(q, S, SA)


# Little stub for testing
if __name__ == "__main__":
	from sys import argv
	filename = argv[1]
	print parse_tpkl(filename)