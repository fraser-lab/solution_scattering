""" This part of the package is for loading data of various types and then 
making traces.

Benjamin Barad
"""
from numpy import recarray, load
from trace import Trace

def parse(filename, mode="tpkl"):
	"""Wrapper function for any loader functions that I may write besides 
	tpkl. Just passes through to the appropriate place based on the `mode` 
	variable."""
	if mode == "tpkl":
		return parse_tpkl(filename)

def parse_tpkl(filename):
	"""Loads tpkl files and generates a corresponding Trace object. Requires
	table.py from the Anfinrud lab, which we will not distribute.
	"""
	try:
		from table import table
	except ImportError:
		print """You do not have the required code accessible to parse tpkl
		files. Try using a different file format for input"""
		raise
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