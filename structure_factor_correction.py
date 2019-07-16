"""
Correct time-resolved data for packing effects (scattering from intermolecular interactions). This script divides all
time resolved DATs by the structure packing factor DAT, and then writes new DAT files with the corrected data.

Usage:
python3    structure_factor_correction.py    sample_directory    structure_packing_factor.dat

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
script, tr_dir, spf_file = argv

### organize files
tr_dir = pathlib.Path(tr_dir)
spf_file = pathlib.Path(spf_file)
tr_diff_files = list(tr_dir.glob(pattern='*diff*.dat'))
avg_off_protein_only_file = str(list(tr_dir.glob(pattern='*protein_only*.dat'))[0])

### read files
spf = pd.read_table(spf_file,names=["q","SA","sigSA"], delim_whitespace=True, skiprows=1)
avg_off = pd.read_table(avg_off_protein_only_file,names=["q","SA","sigSA"], delim_whitespace=True, skiprows=1)

### correct the average-off protein-only curve for packing effects
spf_corrected_int = avg_off.SA/spf.SA
spf_corrected_error = spf_corrected_int*np.sqrt((avg_off.sigSA/avg_off.SA)**2 + (spf.sigSA/spf.SA)**2) ### https://terpconnect.umd.edu/~toh/models/ErrorPropagation.pdf
avg_off_spf_corrected = Trace(avg_off.q, np.empty_like(avg_off.q), np.empty_like(avg_off.q), spf_corrected_error, spf_corrected_int, np.empty_like(avg_off.q))

### scale back up to proper intensities
avg_off_spf_corrected.scale(avg_off, qmin=0.03, qmax=0.1, approach="projection")
avg_off_spf_corrected_scaled = Trace(avg_off.q, np.empty_like(avg_off.q), np.empty_like(avg_off.q), avg_off_spf_corrected.scaled_sigSA, avg_off_spf_corrected.scaled_SA, np.empty_like(avg_off.q))
avg_off_spf_corrected_scaled.write_dat(avg_off_protein_only_file.replace('.dat','spf-corrected.dat'))
# plt.plot(avg_off.q,avg_off.SA, label='avg_off')
# plt.plot(avg_off_spf_corrected.q,avg_off_spf_corrected.SA, label='avg_off_spf_corr')
# plt.plot(avg_off_spf_corrected_scaled.q,avg_off_spf_corrected_scaled.SA, label='avg_off_spf_corr_scaled')

### add corrected protein-only signal to each TR difference signal
for file in tr_diff_files:
	name = file.name
	tr_diff = parse(name)
	tr_sum = tr_diff.add(avg_off_spf_corrected_scaled)
	tr_sum.write_dat(name.replace('.dat','_spf-corrected.dat'))

# plt.xscale('log')
# plt.legend()
# plt.show()


