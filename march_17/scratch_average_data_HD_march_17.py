import csv
from os import listdir
from sys import argv
import pickle as pkl
import subprocess

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd





"""Statics"""
CHI_OUTLIER = 1.5
TEMPS = ["14C"]
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"] # 
# TIMES = ["-10.1us", "100us"]
# TIMES = ["-10.1us", "1us", "10us", "100us", "1ms"]
# TIMES = ["-10.1us"]
# TIMES = ["-10us",  "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms", "1.78ms", "3.16ms", "5.62ms", "10ms"] # 
# MEGAREPS = 2
REPS = range(5,50)
# PREFIX = "CypA-6"
# PREFIX = "CypA-5"

PREFIX = "CypA-WT-1"
PKL_FILENAME = "WT_13C_protein_full_algebraic.pkl"
DATFILE_PREFIX = "WT_13C_buffer"


from parse import parse_tpkl, alg_scale, lin_regress_scale, integration_scale


length = 0
directories = argv[1:]
reference = parse_tpkl("/Users/benjaminbarad/Documents/170302_APS/CypA-NH-Buffer-1/xray_images/CypA-NH-Buffer-1_17_-10us-10.tpkl")
# reference = parse_tpkl("/Volumes/BAB_AGORA/July_Beamline_Trip/Analysis/common/integration/CypA-S99T/CypA-S99T-Buffer-2/xray_images/CypA-S99T-Buffer-2_17_-10us-14_on.tpkl")
# fig,ax = plt.subplots()
# files = []
all_vectors = []
subtracted_vectors = {i: [] for i in TIMES}	
fig1, ax1 = plt.subplots()
for directory in directories:
	# files = listdir(directory)
	# for index, _ in enumerate()
	# print length
	# for megarep in range(MEGAREPS):
		for i in REPS: 
			for temp in TEMPS:
				for index, time in enumerate(TIMES):
					try: 
						# onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
						# onscale = int(onstring.split()[3])
						on = parse_tpkl("{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX, i+1, time))
						ax1.plot(on.q, on.as_vector())
						on_scaled = alg_scale(reference, on)
						# on_scaled = lin_regress_scale(reference, on)
						# on_scaled = on.scale_isosbestic()[0]
						# on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)).as_vector()[80:]
						if index > 0:
							off_count = "-{}".format(index+1)
						else:
							off_count = ""
						# offstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_off beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
						# offscale = int(offstring.split()[3])
						off = parse_tpkl("{0}/{1}_{2}_-10us{3}.tpkl".format(directory, PREFIX, i+1, off_count))
						off_scaled = alg_scale(reference, off)
						# off_scaled = lin_regress_scale(reference, off)
						# off_scaled = off.scale_isosbestic()[0]
						# print "{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX, i+1, time)
						# 
						subtracted = [(on_scaled[0][j] - off_scaled[0][j], np.sqrt(on_scaled[1][j]**2 + off_scaled[1][j]**2)) for j in range(len(on.q))]
						# subtracted = [(on_scaled[0][j], np.sqrt(on_scaled[1][j]**2)) for j in range(len(on.q))]

						# print "appending"
						subtracted_vectors[time].append(subtracted)
					except:
						pass
						print "one or both of the on/off pairs was tossed:"
						print "{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX, i+1, time)
fig1.savefig("unscaled.png")

def chisquared(var, ref):
	nu = len(var)-1.0
	I, sigma = zip(*var) 
	Iref,sigmaref = zip(*ref)
	chi_squared = 1/nu*sum([(I[i]-Iref[i])**2/(sigmaref[i]**2)  for i in range(len(var))])
	# print chi_squared
	return chi_squared

def chi_outliers(vectors, reference_vector):
	list = [chisquared(i, reference_vector) for i in vectors]
	print np.mean(list), np.std(list)
	outlier_list = []
	for i, val in enumerate(list):
		if val > CHI_OUTLIER:
			outlier_list.append(i)
			# print i, val
	return outlier_list

def combine_vectors_outliers(vectors):
	averaged_vector = []
	for i in range(len(vectors[0])):
		value_list = [v[i] for v in vectors]
		# averaged_vector.append(np.mean(value_list))
		#### Ensemble Weighting - mean = sum(value/sigma^2)/sum(1/sigma^2)
		####											std = sqrt(1/sum(1/sigma^2))
		means, stds = zip(*value_list)
		avg_mean = np.mean(means)
		std_mean = np.std(means)
		std_prop = np.sqrt(sum([j**2 for j in stds]))/(len(stds)-1)
		std_tot = np.sqrt(std_mean**2+std_prop**2)
		if i == 0:
			print "Standard deviation of means: ", std_mean
			print "Propogated standard deviation from data points: ", std_prop
			print "Total propogated standard deviation: ", std_tot
		# avg_mean = sum([means[i]/(stds[i]**2) for i,_ in enumerate(means)])/sum([1/stds[i]**2 for i,_ in enumerate(means)])
		# std_tot = np.sqrt(1/sum([1/stds[i]**2 for i,_ in enumerate(means)]))
		# print on.q[i], avg_mean, std_mean, std_prop, std_tot
		averaged_vector.append((avg_mean, std_tot))
	outlier_list = chi_outliers(vectors, averaged_vector)
	if len(outlier_list) > 0:
		new_vectors = [vector for i,vector in enumerate(vectors) if i not in outlier_list]
		print len(new_vectors)
		averaged_vector = combine_vectors_outliers(new_vectors)
	return averaged_vector



q = on.q
fig, ax = plt.subplots()	
plots = []
averaged_vectors = []				
for time in TIMES:
	print "===="
	print time
	vectors = subtracted_vectors[time]
	averaged_vector = combine_vectors_outliers(vectors)
	# outlier_list = chi_outliers(vectors, averaged_vector)
	# final_averaged = chi_outliers()
		# print on.q[i], np.mean(value_list), np.std(value_list)
	plots.append(ax.errorbar(q, zip(*averaged_vector)[0], yerr=zip(*averaged_vector)[1], label=time)[0])
	averaged_vectors.append((time, averaged_vector))
	
# with open("WT_HD_protein_onoff_stdevs.csv", 'wb') as csvfile:
# 	writer = csv.writer(csvfile, delimiter="\t")
# 	writer.writerow(["q"]+[str.join(i[0],"\t","sig", for i in averaged_vectors])
# 	print ["q"]+[i[0] for i in averaged_vectors]
# 	for i in range(len(averaged_vector)):
# 		writer.writerow([on.q[i]]+[[j[1][i][0],j[1][i][1]] for j in averaged_vectors])


def make_pkl(vectors):
	with open(PKL_FILENAME, "wb") as pklfile:
		pkl.dump(averaged_vectors, pklfile)



ax.set_xscale("log", nonposx='clip')
# plt.legend(plots, loc='upper center', bbox_to_anchor=(0.5, 1.25),
#           ncol=3, fancybox=True, shadow=True)
plt.show()
make_pkl(averaged_vectors)
# make_dats(averaged_vectors)