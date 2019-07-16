
"""
Correct static data for packing effects (scattering from intermolecular interactions). This script divides all
static DATs by the structure packing factor DAT, and then writes new DAT files with the corrected data.
Note that you can use the static_pattern_glob (e.g. *protein_only*.dat)
to select from among all *.dat files within the sample directory.

Usage:
python3    spc2.py    sample_directory    static_pattern_glob    structure_packing_factor.dat

Author(s):
Alexander M. Wolff
"""


import pandas as pd
import numpy as np
import pathlib
from sys import argv
from trace import Trace
from parse import parse
# import matplotlib.pyplot as plt

### take command line input
script, static_dir, pattern, spf_file = argv

### organize files
static_dir = pathlib.Path(static_dir)
spf_file = pathlib.Path(spf_file)
static_files = list(static_dir.glob(pattern="*{}*".format(pattern)))
# avg_off_protein_only_file = str(list(tr_dir.glob(pattern='*protein_only*.dat'))[0])

### read files
spf = pd.read_table(spf_file,names=["q","SA","sigSA"], delim_whitespace=True, skiprows=1)


for file in static_files:
	name = file.name
	orig = parse(name)
	### correct the curve for packing effects
	spf_corrected_int = orig.SA/spf.SA
	spf_corrected_error = spf_corrected_int*np.sqrt((orig.sigSA/orig.SA)**2 + (spf.sigSA/spf.SA)**2) ### https://terpconnect.umd.edu/~toh/models/ErrorPropagation.pdf
	static_spf_corrected = Trace(orig.q, np.empty_like(orig.q), np.empty_like(orig.q), spf_corrected_error, spf_corrected_int, np.empty_like(orig.q))
	# ### scale back up to proper intensities
	static_spf_corrected.scale(orig, qmin=0.03, qmax=0.1, approach="projection")
	static_spf_corrected_scaled = Trace(orig.q, np.empty_like(orig.q), np.empty_like(orig.q), static_spf_corrected.scaled_sigSA, static_spf_corrected.scaled_SA, np.empty_like(orig.q))
	static_spf_corrected_scaled.write_dat(name.replace('.dat','spf-corrected.dat'))

# plt.xscale('log')
# plt.legend()
# plt.show()

