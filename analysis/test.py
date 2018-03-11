import parse
# import pathlib
from scipy.stats import chisquare
import scipy.integrate as integrate
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
plt.style.use('ggplot')
from collections import OrderedDict
import numpy as np
from numpy.linalg import svd
from trace import Trace
from new_analysis import real_space_plotter


reference = parse.parse("/Volumes/LaCie/radial_averages/Lysozyme-apo-Buffer-1/xray_images/Lysozyme-apo-Buffer-1_2_3C_1_23.7us.tpkl")


plt.plot(reference.q,reference.SA)
plt.show()