"""
Using tempurature standards and SVD of the water ring, 
calculate the average temperature of the on and off at each 
timepoint in a t-jump series.

Requires Parse and Trace
"""

import csv
import math
from os import listdir
import subprocess
from sys import argv

#import matplotlib
#matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd

from parse import parse, alg_scale

##############################################
# Statics Information
##############################################
#TEMPS = ["5C", "14C", "-5C", "0C", "10C", "13C", "3C", "8C", "18C"]
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", 
	"7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us",  "133us", 
	"178us", "237us", "316us", "422us", "562us", "750us", "1ms"] # 
# TIMES = ["-10.1us", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms"]
STATIC_REPS = range(32)
STATIC_TEMPS = [3, 8, 13, 18, 23, 28] 
# STATIC_TEMPS = [14,21,28]

##############################################
# Time Resolved Information
##############################################
REPS = range(5,40)

TR_DIRECTORIES = ["/mnt/d/T-jump_CypA_best/March2017/Analysis/WAXS/common/integration/CypA/CypA-WT-1/xray_images/", 
	"/mnt/d/T-jump_CypA_best/March2017/Analysis/WAXS/common/integration/CypA/CypA-WT-2/xray_images/", 
	"/mnt/d/T-jump_CypA_best/March2017/Analysis/WAXS/common/integration/CypA/CypA-WT-3/xray_images/", 
	"/mnt/d/T-jump_CypA_best/March2017/Analysis/WAXS/common/integration/CypA/CypA-WT-4/xray_images/", 
	"/mnt/d/T-jump_CypA_best/March2017/Analysis/WAXS/common/integration/CypA/CypA-NH-1/xray_images/"]
TR_PREFIXES = ["CypA-WT-1", "CypA-WT-2", "CypA-WT-3", "CypA-WT-4", "CypA-NH-1"]

assert len(TR_PREFIXES) == len(TR_DIRECTORIES)
STATIC_PREFIX = "CypA-WT-static-1_offPC0T"





length = 0
static_directory = "/mnt/d/T-jump_CypA_best/March2017/Analysis/WAXS/common/integration/CypA/CypA-WT-static-1/xray_images/"


# fig,ax = plt.subplots()
# files = []
vectors = []
subtracted_vectors = []	

##############################################
# Reference dataset is used to scale all data 
##############################################
reference = parse("/mnt/d/T-jump_CypA_best/March2017/Analysis/WAXS/common/integration/CypA/CypA-WT-static-1/xray_images/CypA-WT-static-1_offPC0T13_20.tpkl")


##############################################
# Load and scale all static temperature data
# Goes through each temp and loads them all in turn.
##############################################
for index, temp in enumerate(STATIC_TEMPS):
	for i in STATIC_REPS:
		static_string = "{0}/{1}{2}_{3}.tpkl".format(static_directory, STATIC_PREFIX, temp, i+1)
		static = parse_tpkl(static_string)
		static_scaled = alg_scale(reference, static)[0]
		vectors.append((static_scaled, static_string))



############################################## 
# Load and scale all t-jump data.
# Store the number you load in for each time-resolved dataset because it is variable (outliers removed)
# Iteration: Dataset (usually temp or mutant variant), then repeat iwthin that dataset, then 
# Load and algebraically scale on and off separately for each time in each rep, and store both.
##############################################
lengths = []
for ind,tr_directory in enumerate(TR_DIRECTORIES):
	PREFIX = TR_PREFIXES[ind]
	length = 0
	for i in REPS:
		#for temp in TEMPS:
			for index, time in enumerate(TIMES):
				try: 
					on_string = "{0}/{1}_{2}_{3}.tpkl".format(tr_directory, PREFIX, i+1, time)
					on = parse(on_string)
					on_scaled = alg_scale(reference, on)[0]
					if index > 0:
							off_count = "-{}".format(index+1)
					else:
							off_count = ""
					off_string = "{0}/{1}_{2}_-10us{3}.tpkl".format(tr_directory, PREFIX, i+1, off_count)
					off = parse(off_string)
					off_scaled = alg_scale(reference, off)[0]
					vectors.append((off_scaled, off_string, on_string))
					vectors.append((on_scaled, off_string, on_string))
					length += 2
				except:
					"""On some of our old datasets, we had a different file naming convention,  so we check here"""
					try: 
						on_string = "{0}/{1}_{2}_{3}_on.tpkl".format(tr_directory, PREFIX, i+1, time)
						on = parse(on_string)
						on_scaled = alg_scale(reference, on)[0]
						if index > 0:
								off_count = "-{}".format(index+1)
						else:
								off_count = ""
						off_string = "{0}/{1}_{2}_-10us{3}_on.tpkl".format(tr_directory, PREFIX, i+1, off_count)
						off = parse(off_string)
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
##############################################
## The code below can break if data is included from multiple trips, since 
## different experiment geometries can change the number of Q bins.
##############################################
q = on.q 
length = len(q) 



# matrix = np.matrix([k[0] for k in vectors]).transpose()
for k in vectors[0:3]:
	print k
	print k[0]
	print k[0][1]
##############################################
# SVD with Q Scaling
##############################################
matrix = np.matrix([[q[j]*k[0][j] for j in range(20,length-700)] for k in vectors]).transpose()


u,s,v = svd(matrix, full_matrices=False)

sum_var = np.sum([i**2 for i in s])
print "Variance captured" 
print [(i**2)/sum_var for i in s][0:10]

# print u.shape
# print s.shape
# print v[0]
# print v[1]
# for vector in v[0:7]:
# 	print vector[0]

# print v.slice(0)

##############################################
# Plot right singular vectors to show how they differ per snapshot, and save them out to a csv.
##############################################
fig, ax = plt.subplots()
ax.axvspan(0,29.5,color="0.9", alpha=0.5, linewidth=0)
with open("timepoints.csv", "wb") as csv_output:
	data_struct = v.tolist()[:8]
	wr = csv.writer(csv_output)
	for index in range(len(data_struct[0])):
		wr.writerow([i[index] for i in data_struct])
for vector in v.tolist()[1:4]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	ax.plot([value for value in vector], "-") # , color="#60BD68"
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


##############################################
#Plot left singular vectors to show their shapes, and save them out to a csv.
##############################################
with open("singular_vectors.csv", "wb") as csv_output:
	q_dat = q[20:-700]
	data_struct =u.transpose().tolist()[:8] 
	wr = csv.writer(csv_output)
	for index, q_val in enumerate(q_dat):
		row = [q_val]
		row.extend([i[index] for i in data_struct])
		wr.writerow(row)

for vector in u.transpose().tolist()[1:4]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	x = q[20:-700]
	ax2.plot(x, [value for value in vector], "-") #color="#60BD68"
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
#plt.show()

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


##############################################
# Plot singular values to show variance accounted for by each vector.
##############################################	
fig3, ax3= plt.subplots()
ax3.plot([np.log(i) for i in s][0:10], "-")
fig3.savefig("Singular_values_WT_HD_Unscaled.png")


##############################################
#Separate out the 2nd right singular vector data by temperature and plot as a function of temoperature.
##############################################
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

##############################################
# Fit a quadratic curve to the v vs T data to generate a standard curve for temperature calculation. 
##############################################
m,n,b = np.polyfit(x,y,2, w=y_weight)
print m, n, b

ax4.plot(x, [m*i**2+n*i+b for i in x], ls="--")


##############################################
# Also calculate a linaer fit, but that doesn't fit nearly as well.
##############################################
m1, b1 = np.polyfit(x,y,1)
print m1, b1
ax4.plot(x,[m1*i+b1 for i in x])

fig4.savefig("temperature_fits.png")


##############################################
# Save out standard curve information to a csv
##############################################
with open("temperature_fits.csv", "wb") as csv_output:
	wr = csv.writer(csv_output)
	wr.writerow(["Temp", "Measured", "Quadratic Estimate", "Linear Estimate"])
	for index, x_val in enumerate(x):
		row = [x_val, y[index], m*x_val**2 + n*x_val + b, m1*x_val + b1]
		wr.writerow(row)



##############################################
# Convert time strings to values usable for plotting.
##############################################
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


##############################################
# Iterate through datasets and times and average the ons and offs separately, 
# then calculate temperatures for each based on the SVD quadratic fit.
#
# Plot delta temperatures as a function of time delay for each experiment.
##############################################
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