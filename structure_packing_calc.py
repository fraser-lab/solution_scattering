
import pandas as pd
import numpy as np
# from scipy.stats import linregress
import scipy.optimize
# import matplotlib.pyplot as plt
import pathlib
from sys import argv
from collections import defaultdict
from trace import Trace

script, samp_dir = argv

samp_dir = pathlib.Path(samp_dir)
samp_files = list(samp_dir.glob(pattern='*.dat'))
TEMPS = []
CONCS = []
master  = defaultdict(lambda: defaultdict())

for file in samp_files:
	name = file.name
	parent = file.parent
	samp, data_type, temp, conc, n1, n2, n3 = name.split('_')
	TEMPS.append(temp)
	CONCS.append(conc)
	data = pd.read_table(file,names=["q","SA","sigSA"], delim_whitespace=True, skiprows=1)
	data = data[data.q>=0.03]
	master[temp][conc]=data

TEMPS = sorted(list(set(TEMPS)))
CONCS = sorted(list(set(CONCS)))

def linear(x,a,b):
	return a*x+b

for temp in TEMPS:
	little = [master[temp][CONCS[0]]['q']]
	for conc in CONCS:
		little.append(master[temp][conc]['SA'])

	little.append(master[temp]['PC0']['sigSA'])
	new = pd.concat(little, axis=1, keys=['q','I_pc0', 'I_pc1', 'I_pc2', 'sigSA_pc0'])
	I_0s = []
	A_s = []
	I_0s_errors = []
	for i in new.index:
	    y = np.array([1/(new["I_pc0"][i]),1/(new["I_pc1"][i]),1/(new["I_pc2"][i])])
	    x = np.array([0.050,0.050/3,0.050/9])
	    # model = linregress(x,y)
	    popt, pcov = scipy.optimize.curve_fit(linear, x, y, method='lm', p0=[-10,10], maxfev=50000)

	    # I_0 = np.exp(popt[1])
	    I_0 = 1/popt[1]
	    slope = popt[0]
	    I_0_error = np.sqrt(pcov[1][1])
	    MW = 18500
	    A = slope*I_0/(2*MW)
	    I_0s.append(I_0)
	    A_s.append(A)
	    I_0s_errors.append(I_0_error)

	iz = pd.DataFrame(I_0s)
	iz.index = iz.index + 4
	az = pd.DataFrame(A_s)
	az.index = az.index + 4
	ierrorz = pd.DataFrame(I_0s_errors)
	ierrorz.index = ierrorz.index + 4

	new["I_0"] = iz
	new["A"] = az
	new["I_0s_errors"] = ierrorz
	spf = new.I_pc0/new.I_0
	spf_error = spf*np.sqrt((new.sigSA_pc0/new.I_pc0)**2 + (new.I_0s_errors/new.I_0)**2) ### https://terpconnect.umd.edu/~toh/models/ErrorPropagation.pdf
	structure_packing_factor = Trace(new.q, np.empty_like(new.q), np.empty_like(new.q), spf_error, spf, np.empty_like(new.q))
	structure_packing_factor.write_dat(samp+"_"+temp+"_spf"+".dat")

# print(new)

# plt.figure(num=1, figsize=(9,6))
# plt.scatter(new.q,new.A)
# plt.xlabel("q")
# plt.ylabel("A")
# plt.ylim(-0.001,.001)
# plt.xlim(0,1)

# plt.figure(num=2, figsize=(9,6))
# plt.scatter(new.q,new.I_0)
# # plt.scatter(new.q,new.I_pc0)
# # plt.scatter(new.q,new.I_pc1)
# # plt.scatter(new.q,new.I_pc2)
# plt.xlabel("q")
# plt.ylabel("I_0")
# plt.xscale("log")
# plt.show()
# # plt.ylim(-1,1)