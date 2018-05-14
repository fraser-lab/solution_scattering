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
# from table import table



def parse(filename):
    """Wrapper function for any loader functions that I may write besides 
    tpkl. Just passes through to the appropriate place based on the `mode` 
    variable."""
    if filename.endswith("tpkl"):
        return parse_tpkl_2(filename)
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
    TPKL_HEADER_BYTES = 279 ### this value could vary...original value
    # TPKL_HEADER_BYTES = 290 ### march 2018
    with open(filename, "rb") as f:
        f.seek(TPKL_HEADER_BYTES)
        data = np.fromfile(f, dtype=dt)
        d2 = DataFrame.from_records(data)
    return Trace(d2.q, d2.sigS, d2.S, d2.sigSA, d2.SA, d2.Nj)

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



