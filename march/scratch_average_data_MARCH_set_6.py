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
TEMPS = ["14C", "14.001C"]
TIMES = ["-1us",  "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms", "1.78ms", "3.16ms", "5.62ms", "10ms"] # 
MEGAREPS = 7
REPS = 3
PREFIX = "CypA-6"
# PREFIX = "CypA-5"

from parse import parse_tpkl, alg_scale, Q


length = 0
directories = argv[1:]
# reference = parse_tpkl("/Volumes/BAB_AGORA/CypA-time-series-5/xray_images/CypA-5_2_7.001C_4_100us_off.tpkl".format(directories[0]))
reference = parse_tpkl("/Volumes/BAB_AGORA/CypA-buffer-time-series-7/xray_images/CypA-Buffer-7_4_14.001C_1_56.2us_off.tpkl".format(directories[0]))
# fig,ax = plt.subplots()
# files = []
all_vectors = []
subtracted_vectors = {i: [] for i in TIMES}	
for directory in directories:
	# files = listdir(directory)
	# for index, _ in enumerate()
	# print length
	for megarep in range(MEGAREPS):
		for i in range(REPS):
			for temp in TEMPS:
				for index, time in enumerate(TIMES):
					try: 
						# onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
						# onscale = int(onstring.split()[3])
						on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_{5}_on.tpkl".format(directory, PREFIX, megarep+1, temp, i+1, time))
						on_scaled = alg_scale(reference, on)
						# on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)).as_vector()[80:]
	
						# offstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_off beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
						# offscale = int(offstring.split()[3])
						off = parse_tpkl("{0}/{1}_{2}_{3}_{4}_{5}_off.tpkl".format(directory, PREFIX, megarep+1, temp, i+1, time))
						off_scaled = alg_scale(reference, off)
						print "{0}/{1}_{2}_{3}_{4}_{5}_on.tpkl".format(directory, PREFIX, megarep+1, temp, i+1, time)
	
						subtracted = [on_scaled[j] - off_scaled[j] for j in range(len(on_scaled))]

						# print "appending"
						all_vectors.append(off_scaled)
						all_vectors.append(on_scaled)
						subtracted_vectors[time].append(subtracted)
					except:
						pass
						print "one or both of the on/off pairs was tossed:"
						print "{0}/{1}_{2}_{3}_{4}_{5}_off.tpkl".format(directory, PREFIX, megarep+1, temp, i+1, time)

q = on.q
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
	# print averaged_vector
	plots.append(ax.plot(q, averaged_vector, label=time)[0])
	averaged_vectors.append((time, averaged_vector))
	
with open("set_6_protein_averaged_scaled.csv", 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter="\t")
	writer.writerow(["q"]+[i[0] for i in averaged_vectors])
	print ["q"]+[i[0] for i in averaged_vectors]
	for i in range(len(averaged_vector)):
		writer.writerow([Q[i]]+[j[1][i] for j in averaged_vectors])
	


	
ax.set_xscale("log", nonposx='clip')
ax.legend(plots, loc='upper center', bbox_to_anchor=(0.5, 1.25),
          ncol=3, fancybox=True, shadow=True)
plt.show()