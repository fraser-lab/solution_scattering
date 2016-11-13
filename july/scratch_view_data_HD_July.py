import csv
from sys import argv
import pickle as pkl
import subprocess

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd





"""Statics"""
CHI_OUTLIER = 200
TEMPS = ["14C"]
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"] # 
# TIMES = ["-10.1us", "100us"]
TIMES = ["-10.1us", "1us", "10us", "100us", "1ms"]
# TIMES = ["-10.1us"]
# TIMES = ["-10us",  "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms", "1.78ms", "3.16ms", "5.62ms", "10ms"] # 
# MEGAREPS = 2
REPS = range(1,22)
# PREFIX = "CypA-6"
# PREFIX = "CypA-5"
PREFIX = "CypA-WT-1"
PKL_FILENAME = "S99T_proten_chi.pkl"
DATFILE_PREFIX = "WT_protein"


from parse import parse_tpkl, alg_scale, lin_regress_scale


length = 0
directories = argv[1:]
reference = parse_tpkl("/Users/benjaminbarad/Desktop/July_Trip_Data/CypA-S99T/CypA-S99T-Buffer-2/xray_images/CypA-S99T-Buffer-2_17_-10us-14_on.tpkl")
# reference = parse_tpkl("/Volumes/BAB_AGORA/July_Beamline_Trip/Analysis/common/integration/CypA-S99T/CypA-S99T-Buffer-2/xray_images/CypA-S99T-Buffer-2_17_-10us-14_on.tpkl")
# fig,ax = plt.subplots()
# files = []
all_vectors = []
subtracted_vectors = {i: [] for i in TIMES}	
on_vectors = {i: [] for i in TIMES}	
off_vectors = {i: [] for i in TIMES}	
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
						on = parse_tpkl("{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX, i+1, time))
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
						off = parse_tpkl("{0}/{1}_{2}_-10us{3}_on.tpkl".format(directory, PREFIX, i+1, off_count))
						off_scaled = alg_scale(reference, off)
						# off_scaled = lin_regress_scale(reference, off)
						# off_scaled = off.scale_isosbestic()[0]
						# print "{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX, i+1, time)
						# 
						subtracted = [(on_scaled[0][j] - off_scaled[0][j], np.sqrt(on_scaled[1][j]**2 + off_scaled[1][j]**2)) for j in range(len(on.q))]
						# subtracted = [(on_scaled[0][j], np.sqrt(on_scaled[1][j]**2)) for j in range(len(on.q))]

						# print "appending"
						on_vectors[time].append(zip(*on_scaled))
						off_vectors[time].append(zip(*off_scaled))
						subtracted_vectors[time].append(subtracted)
					except:
						pass
						print "one or both of the on/off pairs was tossed:"
						print "{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX, i+1, time)

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

def combine_vectors_outliers(vectors, on_vectors, off_vectors):
	averaged_vector = []
	averaged_on_vector = []

	for i in range(len(vectors[0])):
		value_list = [v[i] for v in vectors]
		on_value_list = [v[i] for v in on_vectors]
		# off_value_list = [v[i] for v in off_vecrs]
		# averaged_vector.append(np.mean(value_list))
		#### Ensemble Weighting - mean = sum(value/sigma^2)/sum(1/sigma^2)
		####											std = sqrt(1/sum(1/sigma^2))
		means, stds = zip(*value_list)
		on_means, on_stds = zip(*on_value_list)
		# off_means, off_stds = zip(*off_value_list)
		
		avg_mean = np.mean(means)
		std_mean = np.std(means)
		std_prop = np.sqrt(sum([j**2 for j in stds]))/(len(stds)-1)
		std_tot = np.sqrt(std_mean**2+std_prop**2)

		# avg_mean_on = np.mean(on_means)
		# std_mean_on = np.std(on_means)
		# std_prop_on = np.sqrt(sum([j**2 for j in on_stds]))/(len(on_stds)-1)
		# std_tot_on = np.sqrt(std_mean_on**2+std_prop_on**2)

		# avg_mean_off = np.mean(off_means)
		# std_mean_off = np.std(off_means)
		# std_prop_off = np.sqrt(sum([j**2 for j in off_stds]))/(len(off_stds)-1)
		# std_tot_off = np.sqrt(std_mean_off**2+std_prop_off**2)
		# avg_mean = sum([means[i]/(stds[i]**2) for i,_ in enumerate(means)])/sum([1/stds[i]**2 for i,_ in enumerate(means)])
		# std_tot = np.sqrt(1/sum([1/stds[i]**2 for i,_ in enumerate(means)]))
		# print on.q[i], avg_mean, std_mean, std_prop, std_tot
		averaged_vector.append((avg_mean, std_tot))
		# averaged_on_vector.append((avg_mean_on, std_tot_on))
		on_vectors_outlier_free = on_vectors
		off_vectors_outlier_free = off_vectors
		sub_vectors_outlier_free = vectors
		# averaged_off_vector.append((avg_mean_off, std_tot_off))
	outlier_list = chi_outliers(vectors, averaged_vector)
	if len(outlier_list) > 0:
		new_vectors = [vector for i,vector in enumerate(vectors) if i not in outlier_list]
		new_on_vectors = [vector for i,vector in enumerate(on_vectors) if i not in outlier_list]
		new_off_vectors = [vector for i,vector in enumerate(off_vectors) if i not in outlier_list]
		print len(new_vectors)
		print len(new_off_vectors)
		sub_vectors_outlier_free, on_vectors_outlier_free, off_vector_outlier_free = combine_vectors_outliers(new_vectors, new_on_vectors, new_off_vectors)
	return sub_vectors_outlier_free, on_vectors_outlier_free, off_vectors_outlier_free




q = on.q
fig, ax = plt.subplots()	
plots = []
averaged_vectors = []		
averaged_on_vectors = []
free_off_vectors = []	
free_on_vectors = []
free_sub_vectors = []
for time in TIMES:
	print "===="
	print time
	vectors = subtracted_vectors[time]
	sub_vectors_outlier_free, on_vectors_outlier_free, off_vectors_outlier_free = combine_vectors_outliers(vectors, on_vectors[time], off_vectors[time])
	# outlier_list = chi_outliers(vectors, averaged_vector)
	# final_averaged = chi_outliers()
		# print on.q[i], np.mean(value_list), np.std(value_list)
	# plots.append(ax.errorbar(q, zip(*averaged_vector)[0], yerr=zip(*averaged_vector)[1], label=time)[0])
	# averaged_vectors.append((time, averaged_vector))
	free_on_vectors.extend(on_vectors_outlier_free)
	free_off_vectors.extend(off_vectors_outlier_free)
	free_sub_vectors.extend(sub_vectors_outlier_free)

# for vector in free_on_vectors:
# 	ax.plot(q, zip(*vector)[0])
# for vector in free_off_vectors:
# 	ax.plot(q, zip(*vector)[0])
for vector in free_sub_vectors:
	ax.plot(q, zip(*vector)[0])

def combine_off_vectors(free_off_vectors):
	combined_off_vector =[]
	for i in range(len(free_off_vectors[0])):
		value_list = [v[i] for v in free_off_vectors]
		off_means, off_stds = zip(*value_list)
		avg_mean_off = np.mean(off_means)
		std_mean_off = np.std(off_means)
		std_prop_off = np.sqrt(sum([j**2 for j in off_stds]))/(len(off_stds)-1)
		std_tot_off = np.sqrt(std_mean_off**2+std_prop_off**2)
		combined_off_vector.append((avg_mean_off, std_tot_off))
	return combined_off_vector


averaged_off_vector = combine_off_vectors(free_off_vectors)
print averaged_off_vector
# with open("WT_HD_protein_onoff_stdevs.csv", 'wb') as csvfile:
# 	writer = csv.writer(csvfile, delimiter="\t")
# 	writer.writerow(["q"]+[str.join(i[0],"\t","sig", for i in averaged_vectors])
# 	print ["q"]+[i[0] for i in averaged_vectors]
# 	for i in range(len(averaged_vector)):
# 		writer.writerow([on.q[i]]+[[j[1][i][0],j[1][i][1]] for j in averaged_vectors])


def make_pkl(vectors):
	with open(PKL_FILENAME, "wb") as pklfile:
		pkl.dump(averaged_vectors, pklfile)


def make_dats(on_vectors, averaged_off_vector):
	# Offs
	with open("datfiles/{}_off.dat".format(DATFILE_PREFIX), "wb") as dat_file:
			for i in range(len(averaged_off_vector)):
				dat_file.write("{} {} {}\n".format(on.q[i], averaged_off_vector[i][0], averaged_off_vector[i][1]))
	# Ons
	for index,time in enumerate(TIMES):
		vector = on_vectors[index][1]
		with open("datfiles/{}_{}_on.dat".format(DATFILE_PREFIX, time), "wb") as dat_file:
			for i in range(len(vector)):
				dat_file.write("{} {} {}\n".format(on.q[i], vector[i][0], vector[i][1]))
	


ax.set_xscale("log", nonposx='clip')
# plt.legend(plots, loc='upper center', bbox_to_anchor=(0.5, 1.25),
#           ncol=3, fancybox=True, shadow=True)
plt.show()
# make_pkl(averaged_vectors)
# make_dats(averaged_on_vectors, averaged_off_vector)