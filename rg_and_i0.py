# Script to calculate radius-of-gyration (Rg) and extrapolated scattering at 0 angle (I0) using Guinier analysis
#
# Usage: python rg_and_i0.py data_path q_min q_max 
# Argument "data_path" is the directory containing scattering curves
# Arguments "q_min" and "q_max" are the boundaries of the desired Guinier region
# 
# Input scattering curves can be specified by changing the pattern argument in the glob function (see additional note below)
#

import pandas as pd
import numpy as np
import scipy.optimize
import pathlib
from sys import argv
from collections import defaultdict
from trace import Trace
from sys import argv


def linear(x,a,b):
	return b-a*x


script, samp_dir, q_min, q_max = argv
samp_dir_path = pathlib.Path(samp_dir)

### Note : change this pattern to capture different dats within a directory
samp_files = list(samp_dir_path.glob(pattern='*corrected*.dat'))

q_min_squared = float(q_min) ** 2
q_max_squared = float(q_max) ** 2


for file in samp_files:
	name = file.name

	data = pd.read_table(file,names=["q","SA","sigSA"], delim_whitespace=True, skiprows=1)

	x = data.q**2
	y = np.log(data.SA)
	data_mask = np.array(x, dtype=bool)
	data_mask[x>q_max_squared]=False
	data_mask[x<q_min_squared]=False
	x_masked = x[data_mask]
	y_masked = y[data_mask]

	popt, pcov = scipy.optimize.curve_fit(linear, x_masked, y_masked, method='lm', p0=[10,10], maxfev=50000)
	I_0 = np.exp(popt[1])
	slope = popt[0]
	I_0_error = np.sqrt(pcov[1][1])
	I_0_error_scaled = I_0 * I_0_error
	Rg = np.sqrt(3 * slope)
	Rg_error = np.sqrt(pcov[0][0])
	Rg_error_scaled = 0.5 * Rg * ( Rg_error / abs(slope) )

	print(name)
	print("Rg  = {:.2f}    +/-   {:.2f}".format(Rg,Rg_error_scaled))
	print("I_0 = {:.2f} +/- {:.2f}".format(I_0,I_0_error_scaled))
	print('\n')












