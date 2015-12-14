import csv
from os import listdir
import subprocess
from sys import argv

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd



"""Statics"""
TEMPS = ["25C"]
TIMES = ["10ns", "100ns", "1us", "10us"]
REPS = 5
PREFIX = "CypA-11_1"

from parse import parse_tpkl

length = 0
directories = argv[1:]
# fig,ax = plt.subplots()
# files = []
vectors = []
subtracted_vectors = {i: [] for i in TIMES}	
for directory in directories:
	# files = listdir(directory)
	# for index, _ in enumerate()
	# print length
	for i in range(REPS):
		for temp in TEMPS:
			for index, time in enumerate(TIMES):
				try: 
					if i < 5:
						offstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
						offscale = int(offstring.split()[3])
						off = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)).apply_i0_scaling(offscale)
						# off = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)).as_vector()[80:]
				
						onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_off beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
						onscale = int(offstring.split()[3])
						on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, i+1, time)).apply_i0_scaling(onscale)

						# on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, i+1, time)).as_vector()[80:]
						print "{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)
					else:
						index_2 = (index+1) % 4
						if index_2 == 0:
							j = i + 1
						else:
							j = i
						offstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_off beamstop-1.log".format(PREFIX, temp, j+1, TIMES[index_2]), shell=True)
						offscale = int(offstring.split()[3])
						off = parse_tpkl("{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, j+1, TIMES[index_2])).apply_i0_scaling(offscale)						
						# off = parse_tpkl("{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, j+1, TIMES[index_2])).as_vector()[80:]
						
						onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_off beamstop-1.log".format(PREFIX, temp, j+1, TIMES[index_2]), shell=True)
						onscale = int(offstring.split()[3])
						on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, j+1, TIMES[index_2])).apply_i0_scaling(onscale)
						# on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, j+1, TIMES[index_2])).as_vector()[80:]
						print "{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, j+1, TIMES[index_2])
					subtracted = on - off
					# print "appending"
					vectors.append(off)
					vectors.append(on)
					subtracted_vectors[time].append(subtracted)
				except:
					pass
					print "one or both of the on/off pairs was tossed:"
					print "{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, i+1, time)
fig, ax = plt.subplots()	
plots = []
averaged_vectors = []				
for time in TIMES:
	vectors = subtracted_vectors[time]
	averaged_vector = []
	for i in range(len(vectors[0])):
		value_list = [v[i] for v in vectors]
		averaged_vector.append(np.mean(value_list))
	print averaged_vector
	plots.append(ax.plot(averaged_vector, label=time)[0])
	averaged_vectors.append((time, averaged_vector))
	
with open("fast_times_sample_averaged.csv", 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter="\t")
	writer.writerow(["x"]+[i[0] for i in averaged_vectors])
	print ["q"]+[i[0] for i in averaged_vectors]
	for i in range(len(averaged_vector)):
		writer.writerow([(i+6)*0.0025]+[j[1][i] for j in averaged_vectors])
	


	
ax.set_xscale("log", nonposx='clip')
plt.legend(plots)
plt.show()