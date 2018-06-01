
import pandas as pd
import numpy as np
# from scipy.stats import linregress
import scipy.optimize
# import matplotlib.pyplot as plt
import pathlib
from sys import argv
from collections import defaultdict
from trace import Trace

def linear(x,a,b):
	return b-a*x

from sys import argv
from trace import Trace

script, samp_dir = argv

samp_dir_path = pathlib.Path(samp_dir)
samp_files = list(samp_dir_path.glob(pattern='*sum*.dat'))


for file in samp_files:
	name = file.name

	data = pd.read_table(file,names=["q","SA","sigSA"], delim_whitespace=True, skiprows=1)

	x = data.q**2
	y = np.log(data.SA)
	data_mask = np.array(x, dtype=bool)
	data_mask[x>0.008]=False
	x_masked = x[data_mask]
	y_masked = y[data_mask]

	popt, pcov = scipy.optimize.curve_fit(linear, x_masked, y_masked, method='lm', p0=[10,10], maxfev=50000)
	I_0 = np.exp(popt[1])
	slope = popt[0]
	I_0_error = np.sqrt(pcov[1][1])
	I_0_error_scaled = I_0 * I_0_error
	MW = 18500
	Rg = np.sqrt(3 * slope)
	Rg_error = np.sqrt(pcov[0][0])
	Rg_error_scaled = 0.5 * Rg * ((3 * Rg_error) / (3 * abs(slope)))

	print(name)
	print("Rg  = {:.2f}    +/-   {:.2f}".format(Rg,Rg_error_scaled))
	print("I_0 = {:.2f} +/- {:.2f}".format(I_0,I_0_error_scaled))
	print('\n')












