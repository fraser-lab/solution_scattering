""" This part of the package is for loading data of various types and then 
making traces.

Benjamin Barad
"""
from numpy import load
from scipy import stats
from trace import Trace

# Q = [0.0175 + 0.0025 * i for i in range(2125)]
# print Q

def parse(filename, mode="tpkl"):
	"""Wrapper function for any loader functions that I may write besides 
	tpkl. Just passes through to the appropriate place based on the `mode` 
	variable."""
	if mode == "tpkl":
		return parse_tpkl(filename)

def parse_tpkl(filename):
	"""Loads tpkl files and generates a corresponding Trace object."""
	data = load(filename)
	q = data.q
	sigS = data.sigS
	S = data.S
	sigSA = data.sigSA
	SA = data.SA
	Nj = data.Nj
	return Trace(q, sigS, S, sigSA, SA, Nj)

def alg_scale(ref, var):
	SA_ref = ref.SA
	SA_var = var.SA
	# q = SA_ref.q
	# print q
	# return SA_var
	# top = sum([SA_ref[i]*Q[i]*SA_var[i]*Q[i] for i in range(len(SA_ref))])
	top = sum([SA_ref[i]*ref.q[i]*SA_var[i]*ref.q[i] for i in range(len(ref.q))]) # 2074 previously!
	# bottom = sum([(SA_var[i]*Q[i])**2 for i in range(len(SA_var))])
	bottom = sum([(SA_var[i]*ref.q[i])**2 for i in range(len(ref.q))]) # 552,633 # 2074 previously!
	scalar = top/bottom
	print "scalar: ", scalar
	SA_adjusted = [i*scalar for i in SA_var]
	sig_SA_adjusted = [i*scalar for i in var.sigSA]
	return SA_adjusted, sig_SA_adjusted
	# return var.SA, var.sigSA

def lin_regress_scale(ref, var):
	## Performs like algebraic scaling but trains on a much larger q range so alg_scale is preferable.
	slope = stats.linregress([ref.SA[i]*ref.q[i] for i in range(512,593)], [var.SA[i]*ref.q[i] for i in range(512, 593)])[0]
	print slope
	var.SA_adjusted = [i/slope for i in var.SA]
	var.sigSA_adjusted = [i/slope for i in var.sigSA]
	return var.SA_adjusted, var.sigSA_adjusted

# Little stub for testing
if __name__ == "__main__":
	from sys import argv
	filename = argv[1]
	trace = parse_tpkl(filename)
	print trace