""" This part of the package is for loading data of various types and then 
making traces.

Benjamin Barad
"""
from numpy import load
from scipy import stats
from trace import Trace
from pandas import read_table,DataFrame
import numpy as np
import re

# Q = [0.0175 + 0.0025 * i for i in range(2125)]
# print Q

def parse(filename):
    """Wrapper function for any loader functions that I may write besides 
    tpkl. Just passes through to the appropriate place based on the `mode` 
    variable."""
    if filename.endswith("tpkl"):
        return parse_tpkl(filename)
    elif filename.endswith("dat"):
        return parse_dat(filename)
    else:
        raise TypeError('scattering data can only be read from the following filetypes *.tpkl, *.dat')

def parse_dat(filename):
    data = read_table(filename, delimiter="    ", engine='python', skiprows=1, names=['q','I','sigI'])
    q = data.q
    SA = data.I
    sigSA = data.sigI
    S = np.empty_like(data.q)
    sigS = np.empty_like(data.q)
    Nj = np.empty_like(data.q)
    # return q,SA,sigSA
    # sigS = data.sigI
    # S = data.I
    # Nj = data.q
    return Trace(q, sigS, S, sigSA, SA, Nj)


dt = np.dtype({'names': ['q','S','sigS','SA','sigSA','Nj'],
                    'formats': ['<f8','<f8','<f8','<f8','<f8','<i8']})


def parse_tpkl_2(filename):
    """Loads tpkl files and generates a corresponding Trace object.
    """
    TPKL_HEADER_BYTES = 279 ### this value could vary
    with open(filename, "rb") as f:
        f.seek(TPKL_HEADER_BYTES)
        data = np.fromfile(f, dtype=dt)
    return DataFrame.from_records(data)

def parse_tpkl(filename):
    """Loads tpkl files and generates a corresponding Trace object.
    """
    TPKL_HEADER_BYTES = 279 ### this value could vary
    with open(filename, "rb") as f:
        f.seek(TPKL_HEADER_BYTES)
        data = np.fromfile(f, dtype=dt)
    q = data['q']
    sigS = data['sigS']
    S = data['S']
    sigSA = data['sigSA']
    SA = data['SA']
    Nj = data['Nj']
    return Trace(q, sigS, S, sigSA, SA, Nj)


def parse_tpkl_depreciated(filename):
	"""Loads tpkl files and generates a corresponding Trace object.
    This variation is dependent upon table.py from the Anfinrud lab
    because the tpkl files are custom recarray objects"""
	data = load(filename)
	q = data.q
	sigS = data.sigS
	S = data.S
	sigSA = data.sigSA
	SA = data.SA
	Nj = data.Nj
	return Trace(q, sigS, S, sigSA, SA, Nj)


###scaling moved to Trace methods, please use there.

# def alg_scale(ref, var):
#     """Scale by projection of a vector onto a reference vector and determining the magnitude difference."""
#     SA_ref = ref.SA

#     SA_var = var.SA
#     q = ref.q
#     q_SA_ref = SA_ref*q
#     q_SA_var = SA_var*q
#     # q = SA_ref.q
#     # print q
#     # return SA_var
#     # top = sum([SA_ref[i]*Q[i]*SA_var[i]*Q[i] for i in range(len(SA_ref))])
#     top = np.dot(q_SA_var,q_SA_ref) # 2074 previously!
#     # bottom = sum([(SA_var[i]*Q[i])**2 for i in range(len(SA_var))])
#     bottom = np.dot(q_SA_var,q_SA_var) # 552,633 # 2074 previously!
#     scalar = top/bottom
#     print("scalar: {}".format(scalar))
#     SA_adjusted = SA_var * scalar
#     sig_SA_adjusted = var.sigSA * scalar
#     return SA_adjusted, sig_SA_adjusted
#     # return var.SA, var.sigSA

# def integration_scale(ref, var):
#     """Scale by the total number of scattered photons"""
#     SA_ref = ref.SA
#     SA_var = var.SA
#     q = var.q
#     top = np.dot(SA_ref, q)
#     bottom =  np.dot(SA_var, q)
#     scalar = top/bottom
#     SA_adjusted = SA_var * scalar
#     sig_SA_adjusted = var.sigSA * scalar
#     return SA_adjusted, sig_SA_adjusted

def lin_regress_scale(ref, var):
	## Performs like algebraic scaling but trains on a much larger q range so alg_scale is preferable.
	slope = stats.linregress([ref.SA[i]*ref.q[i] for i in range(512,593)], [var.SA[i]*ref.q[i] for i in range(512, 593)])[0]
	print(slope)
	var.SA_adjusted = [i/slope for i in var.SA]
	var.sigSA_adjusted = [i/slope for i in var.sigSA]
	return var.SA_adjusted, var.sigSA_adjusted

def unit_sort(text):
    if text.startswith("-"):
        return 0, text

    elif text.endswith("ns"):
        return 1, text

    elif text.endswith("us"):
        return 2, text

    elif text.endswith("ms"):
        return 3, text

def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]

def tryint(s):
    try:
        return int(s)
    except:
        return s
     
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ] 


def times_numeric(text):
    number = float(text[:-2])
    if text.endswith("ns"):
        return number
    elif text.endswith("us"):
        return 1e3*number
    elif text.endswith("ms"):
        return 1e6*number
    else:
        print("scale could not be calculated")

# Little stub for testing
if __name__ == "__main__":
	from sys import argv
	filename = argv[1]
	trace = parse(filename)
	print(trace)



