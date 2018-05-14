import pandas as pd
import numpy as np
import pathlib
from sys import argv
from trace import Trace

script, tr_dir, spf_file = argv

tr_dir = pathlib.Path(tr_dir)
spf_file = pathlib.Path(spf_file)
tr_sum_files = list(tr_dir.glob(pattern='*sum*.dat'))
tr_diff_files = list(tr_dir.glob(pattern='*diff*.dat'))

spf = pd.read_table(spf_file,names=["q","SA","sigSA"], delim_whitespace=True, skiprows=1)

for file in tr_sum_files:
	name = file.name
	parent = file.parent
	samp, data_type, time = name.split('_')
	time = time.replace('.dat','')

	data = pd.read_table(file,names=["q","SA","sigSA"], delim_whitespace=True, skiprows=1)
	# data = data[data.q>=0.03]
	spf_corrected_int = data.SA/spf.SA
	spf_corrected_error = spf_corrected_int*np.sqrt((data.sigSA/data.SA)**2 + (spf.sigSA/spf.SA)**2) ### https://terpconnect.umd.edu/~toh/models/ErrorPropagation.pdf
	spf_corrected_object = Trace(data.q, np.empty_like(data.q), np.empty_like(data.q), spf_corrected_error, spf_corrected_int, np.empty_like(data.q))
	spf_corrected_object.write_dat(samp+"_"+data_type+"_"+time+"_spf-corrected"+".dat")
