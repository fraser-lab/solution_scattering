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
for directory in directories:
	files = listdir(directory)
	# for index, _ in enumerate()
	# print length
	for file in files:
		vectors.append(parse_tpkl("{0}/{1}".format(directory,file)).as_vector())
		
subtracted_vectors = []	
fig,ax = plt.subplots()
for index, _ in enumerate(vectors[1::2]):
	# subtracted_vectors.append(vectors[2*index]-vectors[2*index-1])
	ax.plot(vectors[2*index]-vectors[2*index-1], "-")
	
plt.show()

# # matrix = np.matrix(vectors).transpose()
# matrix = np.matrix(subtracted_vectors).transpose()

# u,s,v = svd(matrix, full_matrices=False)

# # print u.shape
# # print s.shape
# # print v[0]
# # print v[1]
# # for vector in v[0:7]:
# # 	print vector[0]

# # print v.slice(0)

# fig, ax = plt.subplots()
# i = 0
# for vector in v.tolist()[0:7]:
# 	# print vector
# 	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
# 	ax.plot(range(len(subtracted_vectors)),[value+i*0.3 for value in vector], "-")
# 	i+=1


# # fig2, ax2 = plt.subplots()
# # i = 0
# # for vector in u.transpose().tolist()[0:7]:
# # 	# print vector
# # 	# ax.plot(range(len(vectors)), [value+i for value in vector], "-")
# # 	ax2.plot([value+i for value in vector], "-")
# # 	i+=0.3
	
# # fig3, ax3= plt.subplots()
# # ax3.plot([np.log(i) for i in s][0:10], "-")


# plt.show()
