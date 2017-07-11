import csv
from os import listdir
from sys import argv
import pickle as pkl
import subprocess

# import matplotlib
# matplotlib.use("MacOSX")
# from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd

from time import clock

from itertools import chain

"""Statics"""
CHI_OUTLIER = 1.5
TEMPS = ["20C"]
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms", '3.16ms', '10ms', '31.6ms', '100ms']# 
# TIMES = ["-10.1us", "100us"]
# TIMES = ["-10.1us", "1us", "10us", "100us", "1ms"]
# TIMES = ["-10.1us"]
# TIMES = ["-10us",  "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms", "1.78ms", "3.16ms", "5.62ms", "10ms"] # 
# MEGAREPS = 2
# REPS = range(5,15)
REPS1 = range(5,50)
REPS2 = range(5,50)
# PREFIX = "CypA-6"
# PREFIX = "CypA-5"

### buffer
PREFIX = "Trypsin-APMSF-2"
PKL_FILENAME = "Trypsin-APMSF-2_first-ten_full_algebraic_offs.pkl"
DATFILE_PREFIX = "Trypsin-APMSF-2"

### protein
# PREFIX = "Trypsin-PABA-1"
# PKL_FILENAME = "Trypsin-PABA-1_full_algebraic.pkl"
# DATFILE_PREFIX = "Trypsin-PABA-1"


from parse import parse_tpkl, parse_tpkl_2, alg_scale, lin_regress_scale, integration_scale


length = 0
directories = argv[1:]
reference = parse_tpkl_2("/Volumes/beryllium/Trypsin/Trypsin-BA-Buffer-1/xray_images/Trypsin-BA-Buffer-1_26_-10us-10.tpkl")
# reference = parse_tpkl("/Volumes/BAB_AGORA/July_Beamline_Trip/Analysis/common/integration/CypA-S99T/CypA-S99T-Buffer-2/xray_images/CypA-S99T-Buffer-2_17_-10us-14_on.tpkl")
# fig,ax = plt.subplots()
# files = []
all_vectors = []
subtracted_vectors = {i: [] for i in TIMES} 
# fig1, ax1 = plt.subplots()
for directory in directories:
    # files = listdir(directory)
    # for index, _ in enumerate()
    # print length
    # for megarep in range(MEGAREPS):
    t0_p1 = clock()
    for i in REPS1: 
        for temp in TEMPS:
            for index, time in enumerate(TIMES):
                # try: 
                    # onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
                    # onscale = int(onstring.split()[3])
                    if index > 0:
                        off_count = "-{}".format(index+1)
                    else:
                        off_count = ""

                    on = parse_tpkl("{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX, i+1, time))
                    
                    on_scaled = alg_scale(reference, on)
                    # ax1.plot(on.q, on_scaled[0],c='r')
                    # on_scaled = lin_regress_scale(reference, on)
                    # on_scaled = on.scale_isosbestic()[0]
                    # on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)).as_vector()[80:]
                    
                # except:
                #     pass
                #     print "one or both of the on/off pairs was tossed:"
                #     print "{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX, i+1, time)
    t1_p1 = clock()
    

    t0_p2 = clock()
    for i in REPS2: 
        for temp in TEMPS:
            for index, time in enumerate(TIMES):
                # try:
                    # offstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_off beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
                    # offscale = int(offstring.split()[3])
                    if index > 0:
                        off_count = "-{}".format(index+1)
                    else:
                        off_count = ""

                    off = parse_tpkl_2("{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX, i+1, time))
                    off_scaled = alg_scale(reference, off)

                # except:
                #     pass
                #     print "one or both of the on/off pairs was tossed:"
                #     print "{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX, i+1, time)

    t1_p2 = clock()
    print("FS parser took {} seconds".format(t1_p1-t0_p1))
    print("AW parser took {} seconds".format(t1_p2-t0_p2))
                    
 





