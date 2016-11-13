from os import listdir
from sys import argv

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd



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
print files
subtracted_vectors = []	

for index, _ in enumerate(vectors[1::2]):
	# print files[2*index]
	print files[2*index+1] + "-" + files[2*index] 
	subtracted_vectors.append(vectors[2*index+1]-vectors[2*index])
	# ax.plot(vectors[2*index]-vectors[2*index-1], "-")
# plt.show()

# matrix = np.matrix(vectors).transpose()
matrix = np.matrix(vectors).transpose()

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
for vector in v.tolist()[0:3]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	ax.plot([value+i*0.3 for value in vector], "-")
	i+=1


fig2, ax2 = plt.subplots()
i = 0
for vector in u.transpose().tolist()[0:7]:
	# print vector
	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
	x = [i*0.025 for i in range(len(vector))]	
	ax2.plot(x, [value+i for value in vector], "-")
	i+=0.3
	
fig3, ax3= plt.subplots()
ax3.plot([np.log(i) for i in s][0:10], "-")


plt.show()
