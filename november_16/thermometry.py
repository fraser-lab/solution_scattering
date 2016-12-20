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
TEMPS = ["5C"]
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"] # 
# TIMES = ["-10.1us", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms"]
STATIC_REPS = range(32)
STATIC_TEMPS = [-5,0,5,7,10,12,14,17,22,26] 
# STATIC_TEMPS = [14,21,28]
REPS = range(5,45)
# PREFIX = "CypA-6"
# PREFIX = "CypA-5"
PREFIX = "CypA-Buffer-1"
# PREFIX = "CypA-WT-Buffer-1"
STATIC_PREFIX = "CypA-ConcTemp-1_offBT"
# STATIC_PREFIX = "CypA-WT-static-1_offBT"
from parse import parse_tpkl, alg_scale

length = 0
static_directory = argv[1]
tr_directory = argv[2]
# DIRECTORY_PREFIX = "/Volumes/DatumsDepot/2016/mike/APS_20161110/Analysis/WAXS/common/integration/CypA/"


# fig,ax = plt.subplots()
# files = []
vectors = []
subtracted_vectors = []	

reference = parse_tpkl("/Volumes/DatumsDepot/2016/mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-Buffer-1/xray_images/CypA-Buffer-1_19_-10.1us.tpkl")



for index, temp in enumerate(STATIC_TEMPS):
	for i in STATIC_REPS:
		static_string = "{0}/{1}{2}_{3}.tpkl".format(static_directory, STATIC_PREFIX, temp, i+1)
		static = parse_tpkl(static_string)
		static_scaled = alg_scale(reference, static)[0]
		vectors.append((static_scaled, static_string))


# files = listdir(directory)
# for index, _ in enumerate()
# print length
for i in REPS:
	for temp in TEMPS:
		for index, time in enumerate(TIMES):
			try: 
				# onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
				# onscale = int(onstring.split()[3])
				# on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)).apply_i0_scaling(onscale)[80:]
				on_string = "{0}/{1}_{2}_{3}.tpkl".format(tr_directory, PREFIX, i+1, time)
				on = parse_tpkl(on_string)
				on_scaled = alg_scale(reference, on)[0]
				# on_scaled = on.as_vector()
				if index > 0:
						off_count = "-{}".format(index+1)
				else:
						off_count = ""
				# offstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_off beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
				# offscale = int(offstring.split()[3])
				# off = parse_tpkl("{0}/{1}_{2}_{3}_{4}_off.tpkl".format(directory, PREFIX, temp, i+1, time)).apply_i0_scaling(offscale)[80:]
				off_string = "{0}/{1}_{2}_-10us{3}.tpkl".format(tr_directory, PREFIX, i+1, off_count)
				off = parse_tpkl(off_string)
				off_scaled = alg_scale(reference, off)[0]
				# off_scaled=off.as_vector()
				# print "{0}/{1}_{2}_-10us{3}_on.tpkl".format(directory, PREFIX, i+1, off_count)
				# print "appending"
				vectors.append((off_scaled, off_string, on_string))
				vectors.append((on_scaled, off_string, on_string))
			except:
				pass
				print "one or both of the on/off pairs was tossed:"
				print "{0}/{1}_{2}_{3}.tpkl".format(tr_directory, PREFIX, i+1, time)
		
# print files
q = on.q
length = len(q)
# for index, _ in enumerate(vectors[1::2]):
# 	# print files[2*index]
# 	print files[2*index+1] + "-" + files[2*index] 
# 	subtracted_vectors.append(vectors[2*index+1]-vectors[2*index])
# 	# ax.plot(vectors[2*index]-vectors[2*index-1], "-")
# plt.show()


# matrix = np.matrix([k[0] for k in vectors]).transpose()
matrix = np.matrix([[q[j]*k[0][j] for j in range(length)] for k in vectors]).transpose()
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
ax.axvspan(0,29.5,color="0.9", alpha=0.5, linewidth=0)
for vector in v.tolist()[1:2]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	ax.plot([value for value in vector], "-") # , color="#60BD68"
	i+=1
ax.yaxis.set_ticks_position('left') # this one is optional but I still recommend it...
ax.xaxis.set_ticks_position('bottom')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_xlabel("Image")
ax.set_ylabel("Amplitude (a.u.)")
fig.savefig("timepoints.png")

fig2, ax2 = plt.subplots()
j = 0

for vector in u.transpose().tolist()[1:2]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	x = q
	ax2.plot(x, [value for value in vector], "-") #color="#60BD68"
	j+=0.3
ax2.yaxis.set_ticks_position('left') # this one is optional but I still recommend it...
ax2.xaxis.set_ticks_position('bottom')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
ax2.set_xlabel(r"q ($\AA^{-1}$)")
ax2.set_ylabel(r"q$\cdot$I")
ax2.set_xscale("log", nonposx='clip')
fig2.savefig("Singular_vectors_WT_HD_Unscaled.png")

# zero_vector_values_on = v.tolist()[0][1::2]
# zero_vector_values_off = v.tolist()[0][0::2]
# average_on = np.average(zero_vector_values_on)
# std_on = np.std(zero_vector_values_on)
# average_off =  np.average(zero_vector_values_off)
# std_off = np.std(zero_vector_values_off)
# print average_on, average_off, std_on, std_off
# for index, value in enumerate(v.tolist()[0]):
# 	if index % 2 == 0:
# 		if ((value - average_off) / std_off) > 2.5:
# 			print vectors[index][1]
# 			print vectors[index][2]
# 	elif index % 2 == 1:
# 		if ((value - average_on) / std_on) > 2.5:
# 			print vectors[index][1]
# 			print vectors[index][2]


	
fig3, ax3= plt.subplots()
ax3.plot([np.log(i) for i in s][0:10], "-")
fig3.savefig("Singular_values_WT_HD_Unscaled.png")


fig4,ax4 = plt.subplots()
x = STATIC_TEMPS
y = []
y_weight = []
for index, temp in enumerate(STATIC_TEMPS):
	average_value = np.mean(v.tolist()[1][index*len(STATIC_REPS):index*len(STATIC_REPS)+len(STATIC_REPS)])
	std_value = np.std(v.tolist()[1][index*len(STATIC_REPS):index*len(STATIC_REPS)+len(STATIC_REPS)])
	y.append(average_value)
	y_weight.append(1/std_value)
ax4.plot(x,y)

# plt.show()

m,n,b = np.polyfit(x,y,2, w=y_weight)
print m, n, b

ax4.plot(x, [m*i**2+n*i+b for i in x], ls="--")

m1, b1 = np.polyfit(x,y,1)
print m1, b1

ax4.plot(x,[m1*i+b1 for i in x])
fig4.savefig("temperature_fits.png")
# plt.show()

ons = []
offs = []	
for index,_ in enumerate(TIMES[1:-12]):
	offset = index+1
	total = len(TIMES)
	start = len(STATIC_REPS)*len(STATIC_TEMPS)

	for rep in range(len(REPS)):
		ons.append(v.tolist()[1][start+rep*total*2+offset*2+1])
		offs.append(v.tolist()[1][start+rep*total*2+offset*2])

on_mean = np.mean(ons)
print "checking ons"
print min(np.roots([m,n,b-on_mean]))

off_mean = np.mean(offs)
print "checking offs"
print min(np.roots([m,n,b-off_mean]))


# on_temp = (np.mean(ons)-b)/m
# print np.mean(ons), np.std(ons)
# off_temp = (np.mean(offs)-b)/m
# print np.mean(offs), np.std(offs)
# print on_temp, off_temp


# on_temp = m*np.mean(ons)**2 + n* np.mean(ons) + b
# off_temp = m*np.mean(offs)**2 + n* np.mean(offs) + b

# on_temp = m1*np.mean(ons)+b1
# off_temp = m1*np.mean(offs)+b1

# print on_temp, off_temp
# print on_temp-off_temp

# fig, ax = plt.subplots()
# for index, _ in enumerate(v.tolist()[0]):
	
# 	vector = [i *v.tolist()[0][index] for i in u.transpose().tolist()[0]] + [i *v.tolist()[1][index] for i in u.transpose().tolist()[1]]
# 	ax.plot(vector)
# 	# ax.set_title("index: {0}".format(index))
# fig.savefig("Combination1_2.png")
# plt.show()