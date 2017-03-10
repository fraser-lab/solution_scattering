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
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us",  "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"] # 
# TIMES = ["-10.1us", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms"]
STATIC_REPS = range(32)
STATIC_TEMPS = [3, 8, 13, 18, 23, 28] 
# STATIC_TEMPS = [14,21,28]
REPS = range(5,40)
# PREFIX = "CypA-6"
# PREFIX = "CypA-5"
TR_DIRECTORIES = ["/Volumes/OR_Trail/170302_APS/CypA-NH-Buffer-1/xray_images", "/Volumes/OR_Trail/170302_APS/CypA-WT-Buffer-1/xray_images"] 
TR_PREFIXES = ["CypA-NH-Buffer-1", "CypA-WT-Buffer-1"]
# TR_DIRECTORIES = ["/Volumes/DatumsDepot/2016/Mike/APS_20161110/Analysis/WAXS/common/integration/CypA/CypA-Buffer-1/xray_images"]
# TR_PREFIXES = ["CypA-Buffer-1"]
assert len(TR_PREFIXES) == len(TR_DIRECTORIES)
# PREFIX = "CypA-WT-Buffer-1"
STATIC_PREFIX = "CypA-NH-Buffer-static-1_offBT"
# STATIC_PREFIX = "CypA-WT-static-1_offBT"
from parse import parse_tpkl, alg_scale

length = 0
static_directory = argv[1]
# tr_directories = argv[2:]
# DIRECTORY_PREFIX = "/Volumes/DatumsDepot/2016/mike/APS_20161110/Analysis/WAXS/common/integration/CypA/"


# fig,ax = plt.subplots()
# files = []
vectors = []
subtracted_vectors = []	

reference = parse_tpkl("/Volumes/OR_Trail/170302_APS/CypA-NH-Buffer-1/xray_images/CypA-NH-Buffer-1_17_-10us-10.tpkl")



for index, temp in enumerate(STATIC_TEMPS):
	for i in STATIC_REPS:
		static_string = "{0}/{1}{2}_{3}.tpkl".format(static_directory, STATIC_PREFIX, temp, i+1)
		static = parse_tpkl(static_string)
		static_scaled = alg_scale(reference, static)[0]
		vectors.append((static_scaled, static_string))

lengths = []
for index,_ in enumerate(TR_DIRECTORIES):
	tr_directory = TR_DIRECTORIES[index]
	PREFIX = TR_PREFIXES[index]
	length = 0
	for i in REPS:
		for temp in TEMPS:
			for index, time in enumerate(TIMES):
				try: 
					on_string = "{0}/{1}_{2}_{3}.tpkl".format(tr_directory, PREFIX, i+1, time)
					on = parse_tpkl(on_string)
					on_scaled = alg_scale(reference, on)[0]
					# on_scaled = on.as_vector()
					if index > 0:
							off_count = "-{}".format(index+1)
					else:
							off_count = ""
					off_string = "{0}/{1}_{2}_-10us{3}.tpkl".format(tr_directory, PREFIX, i+1, off_count)
					off = parse_tpkl(off_string)
					off_scaled = alg_scale(reference, off)[0]
					vectors.append((off_scaled, off_string, on_string))
					vectors.append((on_scaled, off_string, on_string))
					length += 2
				except:
					try: 
						on_string = "{0}/{1}_{2}_{3}_on.tpkl".format(tr_directory, PREFIX, i+1, time)
						on = parse_tpkl(on_string)
						on_scaled = alg_scale(reference, on)[0]
						# on_scaled = on.as_vector()
						if index > 0:
								off_count = "-{}".format(index+1)
						else:
								off_count = ""
						off_string = "{0}/{1}_{2}_-10us{3}_on.tpkl".format(tr_directory, PREFIX, i+1, off_count)
						off = parse_tpkl(off_string)
						off_scaled = alg_scale(reference, off)[0]
						vectors.append((off_scaled, off_string, on_string))
						vectors.append((on_scaled, off_string, on_string))
						length += 2
					except:
						pass
						print "one or both of the on/off pairs was tossed:"
						print "{0}/{1}_{2}_{3}.tpkl".format(tr_directory, PREFIX, i+1, time)
	lengths.append(length)
			
print lengths
q = on.q
length = len(q)



# matrix = np.matrix([k[0] for k in vectors]).transpose()
matrix = np.matrix([[q[j]*k[0][j] for j in range(20,length-800)] for k in vectors]).transpose()
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

for index, vector in enumerate(u.transpose().tolist()[1:2]):
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	x = q[20:-800]
	ax2.plot(x, [value/value1 for value, value1 in vector, u.transpose().to_list()[index]], "-") #color="#60BD68"
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
ax4.plot(x,y, ".")


m,n,b = np.polyfit(x,y,2, w=y_weight)
print m, n, b

ax4.plot(x, [m*i**2+n*i+b for i in x], ls="--")

m1, b1 = np.polyfit(x,y,1)
print m1, b1

ax4.plot(x,[m1*i+b1 for i in x])
fig4.savefig("temperature_fits.png")



times_numeric = []
for time in TIMES:
  number = float(time[:-2])
  scale = time[-2:]
  if scale == "ns":
    times_numeric.append(number)
  elif scale == "us":
    times_numeric.append(1000*number)
  elif scale == "ms":
    times_numeric.append(1000*1000*number)
  else:
    print "scale could not be calculated"

fig5,ax5 = plt.subplots()

starting_time = 0
for experiment, _ in enumerate(TR_DIRECTORIES):
	print TR_DIRECTORIES[experiment]

	start = len(STATIC_REPS)*len(STATIC_TEMPS)+sum([i for position, i in enumerate(lengths) if position < experiment])
	total = len(TIMES)
	differences = []
	highs = []
	lows = []
	for index,_ in enumerate(TIMES[starting_time:]):
		ons = []
		offs = []	
		offset = index+starting_time	
		print "====="
		print TIMES[offset]
		for rep in range(len(REPS)):
			ons.append(v.tolist()[1][start+rep*total*2+offset*2+1])
			offs.append(v.tolist()[1][start+rep*total*2+offset*2])
		on_mean = np.mean(ons)
		high = min(np.roots([m,n,b-on_mean]))
		print "Jumped Temperature: %.2f" % high

		off_mean = np.mean(offs)
		low = min(np.roots([m,n,b-off_mean]))
		print "Starting Temperature: %.2f" % low
		highs.append(high)
		lows.append(low)
		differences.append(high-low)
	np.savetxt("{}_highs.txt".format(experiment), highs)
	np.savetxt("{}_lows.txt".format(experiment), lows)
	np.savetxt("{}_differences.txt".format(experiment), differences)
	ax5.plot(times_numeric[1:], differences[1:], "-", label = TR_DIRECTORIES[experiment])

ax5.legend(loc='upper center')
ax5.set_title("")
ax5.set_xlabel("Time (ns)")
ax5.set_ylabel(r"Calculated Temperature Jump ($^\circ$)")
fig5.savefig("cooling.png")
plt.show()





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