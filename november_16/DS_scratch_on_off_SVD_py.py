from os import listdir
from sys import argv
import math
import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd
import pprint as pprint
import scipy.stats.mstats
import pickle as pkl


'''
from parse import parse_tpkl
vectors = []
length = 0
directories = argv[1:]
# fig,ax = plt.subplots()
files = []
for directory in directories:
	files = listdir(directory)
	# for index, _ in enumerate()
	# print length
	for file in files:
		vectors.append(parse_tpkl("{0}/{1}".format(directory,file)).as_vector())
		# ax.plot(parse_tpkl("{0}/{1}".format(directory,file)).as_vector(), "-")
'''
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"]
REPS = range(1,50)
PREFIX = "CypA-1"
#PREFIX = "CypA-Buffer-1"

from parse import parse_tpkl, alg_scale, lin_regress_scale


length = 0
directories = argv[1:]
reference = parse_tpkl("/Users/student/Desktop/CypA-Buffer-1/xray_images/CypA-Buffer-1_19_-10.1us.tpkl")
# fig,ax = plt.subplots()
# files = []
all_vectors = []
subtracted_vectors = {i: [] for i in TIMES}	
for directory in directories:
	# files = listdir(directory)
	# for index, _ in enumerate()
	# print length
	# for megarep in range(MEGAREPS):
		for i in REPS: 
			for index, time in enumerate(TIMES):
				 # try: 
				# onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
				# onscale = int(onstring.split()[3])
				on = parse_tpkl("{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX, i+1, time))
				on_scaled = alg_scale(reference, on)
				all_vectors.append(on_scaled[0][50:])
				#print on_scaled
				# on_scaled = lin_regress_scale(reference, on)
				# on_scaled = on.scale_isosbestic()[0]
				# on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)).as_vector()[80:]
				print "{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX, i+1, time)
				if index > 0:
					off_count = "-{}".format(index+1)
				else:
					off_count = ""
				# offstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_off beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
				# offscale = int(offstring.split()[3])
				off = parse_tpkl("{0}/{1}_{2}_-10us{3}.tpkl".format(directory, PREFIX, i+1, off_count))
				off_scaled = alg_scale(reference, off)
				all_vectors.append(off_scaled[0][50:])
				# off_scaled = lin_regress_scale(reference, off)
				# off_scaled = off.scale_isosbestic()[0]
				# print "{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX, i+1, time)
				# 
				subtracted = [(on_scaled[0][j] - off_scaled[0][j], np.sqrt(on_scaled[1][j]**2 + off_scaled[1][j]**2)) for j in range(len(on.q))]
				# subtracted = [(on_scaled[0][j], np.sqrt(on_scaled[1][j]**2)) for j in range(len(on.q))]

				# print "appending"
				subtracted_vectors[time].append(subtracted)
			 # except:
				# pass
				# print "one or both of the on/off pairs was tossed:"
				# print "{0}/{1}_{2}_{3}_on.tpkl".format(directory, PREFIX, i+1, time)
# pprint.pprint(all_vectors)

matrix = np.matrix(all_vectors).transpose()

u,s,v = svd(matrix, full_matrices=False)
'''
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

def combine_vectors_outliers(vectors):
	averaged_vector = []
	for i in range(len(vectors[0])):
		value_list = [v[i] for v in vectors]
		# averaged_vector.append(np.mean(value_list))
		#### Ensemble Weighting - mean = sum(value/sigma^2)/sum(1/sigma^2)
		####											std = sqrt(1/sum(1/sigma^2))
		means, stds = zip(*value_list)
		avg_mean = np.mean(means)
		std_mean = np.std(means)
		std_prop = np.sqrt(sum([j**2 for j in stds]))/(len(stds)-1)
		std_tot = np.sqrt(std_mean**2+std_prop**2)
		# avg_mean = sum([means[i]/(stds[i]**2) for i,_ in enumerate(means)])/sum([1/stds[i]**2 for i,_ in enumerate(means)])
		# std_tot = np.sqrt(1/sum([1/stds[i]**2 for i,_ in enumerate(means)]))
		# print on.q[i], avg_mean, std_mean, std_prop, std_tot
		averaged_vector.append((avg_mean, std_tot))
	outlier_list = chi_outliers(vectors, averaged_vector)
	if len(outlier_list) > 0:
		new_vectors = [vector for i,vector in enumerate(vectors) if i not in outlier_list]
		print len(new_vectors)
		averaged_vector = combine_vectors_outliers(new_vectors)
	return averaged_vector
# subv = []
# for time in TIMES:
# 	vectors = subtracted_vectors[time]
# 	print vectors
# pprint.pprint(combine_vectors_outliers(subv))

q = on.q
fig, ax = plt.subplots()	
plots = []
averaged_vectors = []				
for time in TIMES:
	print "===="
	print time
	vectors = subtracted_vectors[time]
	averaged_vector = combine_vectors_outliers(vectors)
	# outlier_list = chi_outliers(vectors, averaged_vector)
	# final_averaged = chi_outliers()
		# print on.q[i], np.mean(value_list), np.std(value_list)
	plots.append(ax.errorbar(q, zip(*averaged_vector)[0], yerr=zip(*averaged_vector)[1], label=time)[0])
	averaged_vectors.append((time, averaged_vector))

def make_pkl(vectors):
	with open("average_vectors.pkl", "wb") as pklfile:
		pkl.dump(averaged_vectors, pklfile)
	
# with open("WT_HD_protein_onoff_stdevs.csv", 'wb') as csvfile:
# 	writer = csv.writer(csvfile, delimiter="\t")
# 	writer.writerow(["q"]+[str.join(i[0],"\t","sig", for i in averaged_vectors])
# 	print ["q"]+[i[0] for i in averaged_vectors]
# 	for i in range(len(averaged_vector)):
# 		writer.writerow([on.q[i]]+[[j[1][i][0],j[1][i][1]] for j in averaged_vectors])
ax.set_xscale("log", nonposx='clip')
# plt.legend(plots, loc='upper center', bbox_to_anchor=(0.5, 1.25),
#           ncol=3, fancybox=True, shadow=True)
plt.show()
make_pkl(averaged_vectors)
'''

fig, ax = plt.subplots()
i = 0
for vector in v.tolist()[0:5]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	ax.plot([value+i*0.3 for value in vector], "-")
	i+=1
fig.savefig("timepoints.png")

fig2, ax2 = plt.subplots()
j = 0

for vector in u.transpose().tolist()[0:5]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	x = [(i+7)*0.0025 for i in range(len(vector))]	
	ax2.plot(x, [value for value in vector], "-")
	j+=0.3
fig2.savefig("Singular_vectors_WT_HD_Unscaled.png")

zero_vector_values_on = v.tolist()[0][1::2]
zero_vector_values_off = v.tolist()[0][0::2]
average_on = np.average(zero_vector_values_on)
std_on = np.std(zero_vector_values_on)
average_off =  np.average(zero_vector_values_off)
std_off = np.std(zero_vector_values_off)
# print average_on, average_off, std_on, std_off
for index, value in enumerate(v.tolist()[0]):
	if index % 2 == 0:
		if ((value - average_off) / std_off) > 2.5:
			print all_vectors[index][1]
			print all_vectors[index][2]
	elif index % 2 == 1:
		if math.fabs(((value - average_on) / std_on)) > 2.5:
			print all_vectors[index][1]
			print all_vectors[index][2]
fig3, ax3= plt.subplots()
ax3.plot([np.log(i) for i in s][0:10], "-")
fig3.savefig("Singular_values_WT_HD_Unscaled.png")

plt.show()

# print u.shape
# print s.shape
# print v[0]
# print v[1]
# for vector in v[0:7]:
# 	print vector[0]

# print v.slice(0)
'''
fig, ax = plt.subplots()
i = 0
for vector in v.tolist()[0:3]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	ax.plot([value+i*0.3 for value in vector], "-")
	i+=1


fig2, ax2 = plt.subplots()
i = 0
for vector in u.transpose().tolist()[0:3]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	x = [i+0.0025+0.02 for i in range(len(vector))]	
	#ax2.set_xscale("log")
	ax2.plot(x, [value+i for value in vector], "-")
	i+=0.3
	
fig3, ax3= plt.subplots()
ax3.plot([np.log(i) for i in s][0:10], "-")


plt.show()
'''