import csv
from os import listdir
from sys import argv
import subprocess

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd



"""Statics"""
TEMPS = ["14C"]
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"] # 
# MEGAREPS = 2
REPS = 25
# PREFIX = "CypA-6"
# PREFIX = "CypA-5"
PREFIX = "CypA-S99T-Buffer-2"

from parse import parse_tpkl, alg_scale, Q


length = 0
directories = argv[1:]
reference = parse_tpkl("/Users/benjaminbarad/Desktop/July_Trip_Data/CypA-S99T/CypA-S99T-Buffer-2/xray_images/CypA-S99T-Buffer-2_17_-10us-14_on.tpkl")
# reference = parse_tpkl("/Volumes/BAB_AGORA/July_Beamline_Trip/Analysis/common/integration/CypA-S99T/CypA-S99T-Buffer-2/xray_images/CypA-S99T-Buffer-2_17_-10us-14_on.tpkl")
# fig,ax = plt.subplots()
# files = []
all_vectors = []
subtracted_vectors = {i: [] for i in TIMES}	
for directory in directories:
	# files = listdir(directory)
	# for index, _ in enumerate()
	# print length
	# for megarep in range(MEGAREPS):
		for i in range(REPS): 
			for temp in TEMPS:
				for index, time in enumerate(TIMES):
					# try: 
						# onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
						# onscale = int(onstring.split()[3])
						on = parse_tpkl("{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX, i+1, time))
						on_scaled = alg_scale(reference, on)
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
						# off_scaled = off.scale_isosbestic()[0]
						print "{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX, i+1, time)
	
						subtracted = [on_scaled[j] - off_scaled[j] for j in range(len(on_scaled))]

						# print "appending"
						all_vectors.append(off_scaled)
						all_vectors.append(on_scaled)
						subtracted_vectors[time].append(subtracted)
					# except:
					# 	pass
					# 	print "one or both of the on/off pairs was tossed:"
					# 	print "{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, i+1, time)

q = on.q
print len(q)
fig, ax = plt.subplots()	
plots = []
averaged_vectors = []				
for time in TIMES:
	vectors = subtracted_vectors[time]
	# print vectors
	averaged_vector = []
	for i in range(len(vectors[0])):
		value_list = [v[i] for v in vectors]
		averaged_vector.append(np.mean(value_list))
	print len(averaged_vector)
	plots.append(ax.plot(q, averaged_vector, label=time)[0])
	averaged_vectors.append((time, averaged_vector))
	
with open("S99T_HD_buffer_water_ring.csv", 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter="\t")
	writer.writerow(["q"]+[i[0] for i in averaged_vectors])
	print ["q"]+[i[0] for i in averaged_vectors]
	for i in range(len(averaged_vector)):
		writer.writerow([Q[i]]+[j[1][i] for j in averaged_vectors])
	


	
ax.set_xscale("log", nonposx='clip')
ax.legend(plots, loc='upper center', bbox_to_anchor=(0.5, 1.25),
          ncol=3, fancybox=True, shadow=True)
plt.show()