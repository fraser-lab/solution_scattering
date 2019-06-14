"""
Parse provides a library of pre-built parse functions for different formats
of solution scattering data. Functions are included which can convert data
stored within the filename into python objects for downstream processing.

Author(s):
Benjamin A. Barad
Alexander M. Wolff
"""
from numpy import load
from scipy import stats
from trace import Trace
from pandas import read_table,DataFrame
import numpy as np
import re
# from table import table



def parse(filename):
    """
    A wrapper function which uses the suffix of the filename
    to determine which method to call upon for parsing the file.
    
    Parameters:
    filename (str): path of file to be analyzed
    
    Returns:
    Trace:custom object built to hold a single scattering curve and
    associated values
    """
    if filename.endswith("tpkl"):
        return parse_tpkl_2(filename)
    elif filename.endswith("dat"):
        return parse_dat(filename)
    else:
        raise TypeError('scattering data can only be read from the following filetypes *.tpkl, *.dat')

def parse_dat(filename):
    """
    A function to parse flat text files, with columns separated by spaces.
    
    Parameters:
    filename (str): path of file to be analyzed
    
    Returns:
    Trace:custom object built to hold a single scattering curve and
    associated values
    """
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
    """
    A function to parse custom recarray objects.
    
    Parameters:
    filename (str): path of file to be analyzed
    
    Returns:
    Trace:custom object built to hold a single scattering curve and
    associated values
    """
    TPKL_HEADER_BYTES = 279 ### this value could vary...original value
    # TPKL_HEADER_BYTES = 290 ### march 2018
    with open(filename, "rb") as f:
        f.seek(TPKL_HEADER_BYTES)
        data = np.fromfile(f, dtype=dt)
        d2 = DataFrame.from_records(data)
    return Trace(d2.q, d2.sigS, d2.S, d2.sigSA, d2.SA, d2.Nj)

def parse_tpkl(filename):
    """
    A function to parse custom recarray objects.
    
    Parameters:
    filename (str): path of file to be analyzed
    
    Returns:
    Trace:custom object built to hold a single scattering curve and
    associated values
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
    """
    A function to parse custom recarray objects.
    This variation is dependent upon table.py from the Anfinrud Lab
    
    Parameters:
    filename (str): path of file to be analyzed
    
    Returns:
    Trace:custom object built to hold a single scattering curve and
    associated values
    """
    data = load(filename)
    q = data.q
    sigS = data.sigS
    S = data.S
    sigSA = data.sigSA
    SA = data.SA
    Nj = data.Nj
    return Trace(q, sigS, S, sigSA, SA, Nj)


def unit_sort(text):
    """
    A function to sort files when timepoints are encoded within the filename
    using common abbreviations.
    
    Parameters:
    text (str): filename e.g. (proteinX_100ns)
    
    Returns:
    Tuple(int, text): for sorting files
    """
    if text.startswith("-"):
        return 0, text

    elif text.endswith("ns"):
        return 1, text

    elif text.endswith("us"):
        return 2, text

    elif text.endswith("ms"):
        return 3, text

def atof(text):
    """
    Function to test whether a str can be rendered as a float.
    Useful when retreiving experimental parameters that were
    stored within the filename.
    
    Parameters:
    text (str)
    
    Returns:
    float or text
    """
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text):
    '''
    A regex function to help with sorting strings containing numeric and non-numeric chunks.
    
    Parameters:
    text (str): e.g. if iterated over this list of strings [a1, a11, a2, a13]
    
    Returns:
    int: e.g. would return a list with ints assigning sort position [1, 3, 2, 4]
    '''
    return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]

def tryint(s):
    """
    Function to test whether a str can be rendered as a int.
    Useful when retreiving experimental parameters that were
    stored within the filename.
    
    Parameters:
    text (str)
    
    Returns:
    int or text
    """
    try:
        return int(s)
    except:
        return s
     
def alphanum_key(s):
    """ 
    Turn a string into a list of string and number chunks.
    "z23a" -> ["z", 23, "a"]
    
    Parameters:
    text (str)
    
    Returns:
    list of str and int values
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ] 


def times_numeric(text):
    """
    A function to convert timepoints encoded within the filename
    into a float corresponding to the value in nanoseconds.
    
    Parameters:
    text (str): e.g. (100us)
    
    Returns:
    float: e.g. (100,000)
    """
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



