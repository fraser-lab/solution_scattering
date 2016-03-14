from os import listdir
from sys import argv

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
# from numpy.linalg import svd

from parse import parse_tpkl


vectors_to_show = ["_data/CypA-11/xray_images/CypA-11_1_25C_6_100ns_on.tpkl", "_data/CypA-11/xray_images/CypA-11_1_25C_6_100ns_off.tpkl","_data/CypA-11/xray_images/CypA-11_1_25C_6_10ns_on.tpkl", "_data/CypA-11/xray_images/CypA-11_1_25C_6_10ns_off.tpkl", "_data/CypA-11/xray_images/CypA-11_1_25C_5_10us_on.tpkl", "_data/CypA-11/xray_images/CypA-11_1_25C_5_10us_off.tpkl"]


fig, ax = plt.subplots()
plots = []
# for vector in vectors_to_show:
# 	print vector
# 	y = parse_tpkl(vector).as_vector()
# 	plots.append(ax.plot(y, label=vector)[0])
	
	
plots.append(ax.plot(parse_tpkl(vectors_to_show[1]).as_vector() - parse_tpkl(vectors_to_show[0]).as_vector(), label=vectors_to_show[1])[0])
plots.append(ax.plot(parse_tpkl(vectors_to_show[3]).as_vector() - parse_tpkl(vectors_to_show[2]).as_vector(), label=vectors_to_show[3])[0])
plots.append(ax.plot(parse_tpkl(vectors_to_show[4]).as_vector() - parse_tpkl(vectors_to_show[5]).as_vector(), label=vectors_to_show[5])[0])
ax.legend(plots)
plt.show()
