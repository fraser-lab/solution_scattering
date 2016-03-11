from os import listdir
from sys import argv
import subprocess

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd

"""Statics"""
TEMPS = ["7C", "7.001C"]
TIMES = ["-1us",  "10ns", "31.6ns", "100ns", "316ns", "1us", "3.16us", "10us", "31.6us", "100us", "316us", "1ms", "3.16ms", "10ms"] # 
REPS = 66
PREFIX = "CypA-3"

from parse import parse_tpkl

length = 0
directories = argv[1:]
# fig,ax = plt.subplots()
# files = []
vectors = []
subtracted_vectors = []	
for directory in directories:
	# files = listdir(directory)
	# for index, _ in enumerate()
	# print length
	for i in range(REPS):
		# for temp in TEMPS:
			for time in TIMES:
					# onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
					# onscale = int(onstring.split()[3])
					# on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)).apply_i0_scaling(onscale)[80:]
					on = parse_tpkl("{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX, i+1, time)).as_vector()

					# offstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_off beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
					# offscale = int(offstring.split()[3])
					# off = parse_tpkl("{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, i+1, time)).apply_i0_scaling(offscale)[80:]
					off = parse_tpkl("{0}/{1}_{2}_{3}_off.tpkl".format(directory, PREFIX, i+1, time)).as_vector()


					print "{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX, i+1, time)

					subtracted = on - off

					# print "appending"
					vectors.append(off)
					vectors.append(on)
					subtracted_vectors.append(subtracted)

					# print "one or both of the on/off pairs was tossed:"
					# print "{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, i+1, time)
			
# print files


# for index, _ in enumerate(vectors[1::2]):
# 	# print files[2*index]
# 	print files[2*index+1] + "-" + files[2*index] 
# 	subtracted_vectors.append(vectors[2*index+1]-vectors[2*index])
# 	# ax.plot(vectors[2*index]-vectors[2*index-1], "-")
# plt.show()

matrix = np.matrix(vectors).transpose()
# matrix = np.matrix(subtracted_vectors).transpose()

u,s,v = svd(matrix, full_matrices=False)

# print u.shape
# print s.shape
# print v[0]
# print v[1]
# for vector in v[0:7]:
# 	print vector[0]

# print v.slice(0)

fig, ax = plt.subplots()
i = 0
for vector in v.tolist()[2:4]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	ax.plot([value for value in vector], "-")

fig.savefig("timepoints.png")

fig2, ax2 = plt.subplots()
j = 0

for vector in u.transpose().tolist()[2:4]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	x = [(i+7)*0.0025 for i in range(len(vector))]	
	ax2.plot(x, [value for value in vector], "-")
	j+=0.3
ax2.set_xscale("log", nonposx='clip')
fig2.savefig("Singular_vectors.png")
	
fig3, ax3= plt.subplots()
ax3.plot([np.log(i) for i in s][0:100], "-")
fig3.savefig("Singular_values.png")

plt.show()

# fig, ax = plt.subplots()
# for index, _ in enumerate(v.tolist()[0]):
	
# 	vector = [i *v.tolist()[0][index] for i in u.transpose().tolist()[0]] + [i *v.tolist()[1][index] for i in u.transpose().tolist()[1]]
# 	ax.plot(vector)
# 	# ax.set_title("index: {0}".format(index))
# fig.savefig("Combination1_2.png")
# plt.show()