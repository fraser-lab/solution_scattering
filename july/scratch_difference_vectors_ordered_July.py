import math
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
# TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"] # 
TIMES = ["-10us",  "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms", "1.78ms", "3.16ms", "5.62ms", "10ms"] # 

MEGAREPS = 2
REPS = 21
# PREFIX = "CypA-6"
# PREFIX = "CypA-5"
PREFIX = "CypA-S99T-2"

from parse import parse_tpkl, alg_scale

length = 0
directories = argv[1:]
# fig,ax = plt.subplots()
# files = []
vectors = []
subtracted_vectors = []	
reference = parse_tpkl("/Users/benjaminbarad/Desktop/July_Trip_Data/CypA-S99T/CypA-S99T-Buffer-2/xray_images/CypA-S99T-Buffer-2_17_-10us-14_on.tpkl")

for directory in directories:
	# files = listdir(directory)
	# for index, _ in enumerate()
	# print length
	for megarep in range(MEGAREPS):
		for i in range(REPS):
			for temp in TEMPS:
				for index, time in enumerate(TIMES):
					if i == 0 and megarep == 0:
						continue
					try: 
						# onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
						# onscale = int(onstring.split()[3])
						# on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)).apply_i0_scaling(onscale)[80:]

						on_string = "{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, megarep+1, i+1, time)
						on = parse_tpkl(on_string)
						# on_scaled = alg_scale(reference, on)
						on_scaled = on.as_vector()
						if index > 0:
								off_count = "-{}".format(index+1)
						else:
								off_count = ""
						# offstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_off beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
						# offscale = int(offstring.split()[3])
						# off = parse_tpkl("{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, i+1, time)).apply_i0_scaling(offscale)[80:]
						# off_string = "{0}/{1}_{2}_-10us{3}_on.tpkl".format(directory, PREFIX, i+1, off_count)
						off_string = "{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, megarep+1, i+1, time)
						off = parse_tpkl(off_string)
						# off_scaled = alg_scale(reference, off)
						off_scaled=off.as_vector()
						
						# print "{0}/{1}_{2}_-10us{3}_on.tpkl".format(directory, PREFIX, i+1, off_count)

						subtracted = [on_scaled[j] - off_scaled[j] for j in range(len(on_scaled))]

						# print "appending"
						vectors.append((off_scaled, off_string, on_string))
						vectors.append((on_scaled, off_string, on_string))
						subtracted_vectors.append(subtracted)
					except:
						pass
						print "one or both of the on/off pairs was tossed:"
						print "{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, i+1, time)
				
# print files


# for index, _ in enumerate(vectors[1::2]):
# 	# print files[2*index]
# 	print files[2*index+1] + "-" + files[2*index] 
# 	subtracted_vectors.append(vectors[2*index+1]-vectors[2*index])
# 	# ax.plot(vectors[2*index]-vectors[2*index-1], "-")
# plt.show()

matrix = np.matrix([i[0] for i in vectors]).transpose()
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
for vector in v.tolist()[0:5]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	ax.plot([value+i*0.3 for value in vector], "-")
	i+=1
fig.savefig("timepoints.png")

# fig2, ax2 = plt.subplots()
# j = 0

# for vector in u.transpose().tolist()[0:5]:
# 	# print vector
# 	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
# 	x = [(i+7)*0.0025 for i in range(len(vector))]	
# 	ax2.plot(x, [value for value in vector], "-")
# 	j+=0.3
# fig2.savefig("Singular_vectors_WT_HD_Unscaled.png")

zero_vector_values_on = v.tolist()[0][1::2]
zero_vector_values_off = v.tolist()[0][0::2]
average_on = np.average(zero_vector_values_on)
std_on = np.std(zero_vector_values_on)
average_off =  np.average(zero_vector_values_off)
std_off = np.std(zero_vector_values_off)
print average_on, average_off, std_on, std_off
for index, value in enumerate(v.tolist()[0]):
	if index % 2 == 0:
		if math.abs((value - average_off) / std_off) > 2.5:
			print vectors[index][1]
			print vectors[index][2]
	elif index % 2 == 1:
		if math.abs((value - average_on) / std_on) > 2.5:
			print vectors[index][1]
			print vectors[index][2]


	
# fig3, ax3= plt.subplots()
# ax3.plot([np.log(i) for i in s][0:10], "-")
# fig3.savefig("Singular_values_WT_HD_Unscaled.png")

plt.show()

# fig, ax = plt.subplots()
# for index, _ in enumerate(v.tolist()[0]):
	
# 	vector = [i *v.tolist()[0][index] for i in u.transpose().tolist()[0]] + [i *v.tolist()[1][index] for i in u.transpose().tolist()[1]]
# 	ax.plot(vector)
# 	# ax.set_title("index: {0}".format(index))
# fig.savefig("Combination1_2.png")
# plt.show()