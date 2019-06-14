from collections import namedtuple
from os import listdir, makedirs, path
import os
import shutil
from sys import argv
import math
import pickle as pkl
import pprint as pprint
import warnings

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd
import scipy.stats.mstats

import parse
######################################
# Constants that may relate to geometry or dataset issues.
######################################
QMIN = 0.03
QMAX = 4.28
SKIP_INDEX = 30
OUTLIER_CUTOFF = 2.5

def make_dir(f):
    if not os.path.exists(f):
        os.makedirs(f)


######################################
# Dataset information for loading in data
######################################

TIMES = [ "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", 
	"7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", 
	"133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"]
REPS = range(0,50)

PREFIX_LIST = ["CypA-WT-1",
				 "CypA-WT-Buffer-1",
				 "CypA-S99T-2",
				 "CypA-S99T-Buffer-2",
				 "CypA-1",
				 "CypA-Buffer-1",
				 "CypA-2", 
				 "CypA-Buffer-2",
				 "CypA-3", 
				 "CypA-Buffer-3",
				 "CypA-4", 
				 "CypA-Buffer-4",
				 "CypA-5", 
				 "CypA-Buffer-5",
				 "CypA-WT-1",
				 "CypA-WT-Buffer-1",
				 "CypA-WT-2",
				 "CypA-WT-Buffer-2",
				 "CypA-WT-3",
				 "CypA-WT-Buffer-3",
				 "CypA-WT-4",
				 "CypA-WT-Buffer-4",
				 "CypA-NH-1",
				 "CypA-NH-Buffer-1",
				 "CypA-H-1",
				 "CypA-N-1"]

directories = ["/Volumes/DatumsDepot/2016/Mike/APS_20160701/July_Beamline_Trip/Analysis/common/integration/CypA-WT/CypA-WT-1/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20160701/July_Beamline_Trip/Analysis/common/integration/CypA-WT/CypA-WT-Buffer-1/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20160701/July_Beamline_Trip/Analysis/common/integration/CypA-S99T/CypA-S99T-2/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20160701/July_Beamline_Trip/Analysis/common/integration/CypA-S99T/CypA-S99T-Buffer-2/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-1/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-Buffer-1/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-2/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-Buffer-2/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-3/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-Buffer-3/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-4/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-Buffer-4/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-5/xray_images",
			 "/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-Buffer-5/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-1/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-Buffer-1/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-2/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-Buffer-2/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-3/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-Buffer-3/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-4/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-Buffer-4/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-NH-1/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-NH-Buffer-1/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-H-1/xray_images",
			 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-N-1/xray_images"]

# PREFIX_LIST = [
# 							 "CypA-WT-1",
#  							 "CypA-WT-Buffer-1",
#  							 "CypA-WT-3"]
# #"/Volumes/DatumsDepot/2016/Mike/APS_20160701/July_Beamline_Trip/Analysis/common/integration/CypA-WT/CypA-WT-1/xray_images",
# directories = [
# 							 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-1/xray_images",
#  							 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-Buffer-1/xray_images",
#  							 "/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-3/xray_images"]


######################################
# Load reference
######################################

length = 0
reference = parse.parse("/Volumes/DatumsDepot/2017/Mike/APS_20170302/Analysis/WAXS/common/integration/CypA/CypA-WT-Buffer-2/xray_images/CypA-WT-Buffer-2_9_-10us-22.tpkl")

# reference = parse.parse("/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-Buffer-4/xray_images/CypA-Buffer-4_9_-10us-22.tpkl")
# fig,ax = plt.subplots()
# files = []
Dataset = namedtuple('Dataset', ['directory', 'onavg', 'onstd', 'offavg', 'offstd'])
datasets = []

## Iterate the whole process over each directory.
for directory_index, directory in enumerate(directories):
	all_vectors = []
	subtracted_vectors = {i: [] for i in TIMES}	
	print "================="
	print "Starting directory {}".format(directory)
	print "================="
	######################################
	# Replace files from old quarantines before 
	######################################
	print "Moving quarantined data out of quarantine dir"
	quarantine_dir = "{}/_QUARANTINE".format(directory)
	print quarantine_dir
	if os.path.exists(quarantine_dir):
		files = os.listdir(quarantine_dir)
		for f in files:
			src = "{}/{}".format(quarantine_dir, f)
			dst = "{}/{}".format(directory, f)
			shutil.move(src,dst)

	print "Copying old quarantined data out of quarantine"
	old_quarantine_dir = "{}/_QUARANTINED_DATA".format(directory)
	print old_quarantine_dir
	if os.path.exists(old_quarantine_dir):
		files = os.listdir(old_quarantine_dir)
		for f in files:
			src = "{}/{}".format(old_quarantine_dir, f)
			dst = "{}/{}".format(directory, f)
			shutil.copy(src,dst)


	######################################
	# Parse, Load, Scale Data
	######################################
	print "Loading Data"
	for i in REPS: 
		for index, time in enumerate(TIMES):
			try:
				onstring = "{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX_LIST[directory_index], i+1, time)
				on = parse.parse(onstring)
				# on_scaled = on.scale(reference, qmin=QMIN, qmax=QMAX, approach="algebraic")				
				on_scaled = on.SA, on.sigSA
				if index > 0:
					off_count = "-{}".format(index+1)
				else:
					off_count = ""

				offstring = "{0}/{1}_{2}_-10us{3}.tpkl".format(directory, PREFIX_LIST[directory_index], i+1, off_count)
				off = parse.parse(offstring)
				# off_scaled = off.scale(reference, qmin=QMIN, qmax=QMAX, approach="algebraic")
				off_scaled = off.SA, off.sigSA
				subtracted = [(on_scaled[0][j] - off_scaled[0][j], np.sqrt(on_scaled[1][j]**2 + off_scaled[1][j]**2)) for j in range(len(on.q))]

				all_vectors.append((off_scaled[0][SKIP_INDEX:], onstring, offstring))
				all_vectors.append((on_scaled[0][SKIP_INDEX:], onstring, offstring))
				subtracted_vectors[time].append(subtracted)
			except:
				try:
					onstring = "{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX_LIST[directory_index], i+1, time)
					on = parse.parse(onstring)
					# on_scaled = on.scale(reference, qmin=QMIN, qmax=QMAX, approach="algebraic")				
					on_scaled = on.SA, on.sigSA
					if index > 0:
						off_count = "-{}".format(index+1)
					else:
						off_count = ""

					offstring = "{0}/{1}_{2}_-10us{3}_on.tpkl".format(directory, PREFIX_LIST[directory_index], i+1, off_count)
					off = parse.parse(offstring)
					# off_scaled = off.scale(reference, qmin=QMIN, qmax=QMAX, approach="algebraic")
					off_scaled = off.SA, off.sigSA
					subtracted = [(on_scaled[0][j] - off_scaled[0][j], np.sqrt(on_scaled[1][j]**2 + off_scaled[1][j]**2)) for j in range(len(on.q))]

					all_vectors.append((off_scaled[0][SKIP_INDEX:], onstring, offstring))
					all_vectors.append((on_scaled[0][SKIP_INDEX:], onstring, offstring))
					subtracted_vectors[time].append(subtracted)
				except:
					print "One or both of the on/off pairs was missing"
					print onstring

	######################################
	# SVD and Outlier Detection
	######################################
	print "Running SVD"
	matrix = np.matrix([i[0] for i in all_vectors]).transpose()
	u,s,v = svd(matrix, full_matrices=False)
	zero_vector_values_on = v.tolist()[0][1::2]
	zero_vector_values_off = v.tolist()[0][0::2]
	average_on = np.average(zero_vector_values_on)
	std_on = np.std(zero_vector_values_on)
	average_off =  np.average(zero_vector_values_off)
	std_off = np.std(zero_vector_values_off)
	print average_on, average_off, std_on, std_off
	datasets.append(Dataset(directory, average_on, std_on, average_off, std_off))
	quarantine_set = set()
	print "Finding Outliers"
	for index, value in enumerate(v.tolist()[0]):
		# print all_vectors[index][1]
		if index % 2 == 0:
			if ((value - average_off) / std_off) > OUTLIER_CUTOFF:
				# print all_vectors[index][1]
				print all_vectors[index][2]
				quarantine_set.add(all_vectors[index][1])
				quarantine_set.add(all_vectors[index][2])
		elif index % 2 == 1:
			if math.fabs(((value - average_on) / std_on)) > OUTLIER_CUTOFF:
				print all_vectors[index][1]
				# print all_vectors[index][2]
				quarantine_set.add(all_vectors[index][1])
				quarantine_set.add(all_vectors[index][2])
	######################################
	# Quarantine Outliers
	######################################
	print "Moving outliers into {}".format(quarantine_dir)
	make_dir(quarantine_dir) # I made this directory above
	if quarantine_set:
		for file in quarantine_set:
			basename =  path.basename(file)
			print basename
			dst = "{}/{}".format(quarantine_dir, basename)
			try:
				shutil.move(file, dst)
			except:
				warnings.warn("Move from {} to {} failed - please revisit manually".format(file,dst))
	quarantine_set.clear()
	del quarantine_set
	######################################
	# Plots
	######################################
	print "Plotting Output"
	fig, ax = plt.subplots()
	i = 0
	for vector in v.tolist()[0:1]:
		# print vector
		# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
		ax.plot([value+i*0.3 for value in vector], "-")
		i+=1
	fig.savefig("timepoints_{}.png".format(PREFIX_LIST[directory_index]))
	fig2, ax2 = plt.subplots()
	j = 0

	for vector in u.transpose().tolist()[0:1]:
		# print vector
		# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
		x = [(i+7)*0.0025 for i in range(len(vector))]	
		ax2.plot(x, [value for value in vector], "-")
		j+=0.3
	fig2.savefig("Singular_vectors_WT_HD_Unscaled_{}.png".format(PREFIX_LIST[directory_index]))


	fig3, ax3= plt.subplots()
	ax3.plot([np.log(i) for i in s][0:10], "-")
	fig3.savefig("Singular_values_WT_HD_Unscaled_{}.png".format(PREFIX_LIST[directory_index]))

with open("dataset_statistics.csv", 'wb') as f:
	f.write("Directory,On Average,On STD,Off Average,Off STD\n")
	for d in datasets:
		f.write("{},{},{},{},{}\n".format(d.directory, d.onavg, d.onstd, d.offavg, d.offstd))




