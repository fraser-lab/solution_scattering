from os import listdir
from sys import argv

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd

from parse import parse

"""Statics"""
TEMP = "25C"
TIMES = ["10ns", "100ns", "1us", "10us"]
BUFFER_DIRECTORY = "_data/CypA-Buffer-12/xray_images/"
SAMPLE_DIRECTORY = "_data/CypA-11/xray_images"

def subtractify(buffer_list, sample_list):
	
	


buffer_dict = {i: [k for k in listdir(BUFFER_DIRECTORY) if (TEMP in k and i in k)] for i in TIMES}
sample_dict = {i: [k for k in listdir(SAMPLE_DIRECTORY) if (TEMP in k and i in k)] for i in TIMES}


data = []
fig, ax = plt.subplots()
for time in TIMES:
	averaged_data = subtractify(buffer_dict[i], sample_dict[i])
	

