"""A script to average together traces with error 
propogation and generate .dat files suitable for 
ATSAS analysis

Benjamin Barad, 7 Apr 2016
"""
import math
from os import listdir
from sys import argv
import subprocess

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd

from parse import parse

### Standard deviation scales linearly with isosbestic scaling: http://math.stackexchange.com/questions/682842/scaling-the-normal-distribution
### Variance goes linearly with averaging independent measurements - not stdev

### 
TEMPS = ["14C"]
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"] # 
# MEGAREPS = 2
REPS = 25
# PREFIX = "CypA-6"
# PREFIX = "CypA-5"
PREFIX = "CypA-WT-Buffer-1"

rootname = argv[1]

ons = {i: [] for i in TIMES}
offs = {i: [] for i in TIMES}

# files = listdir(directory)
# for index, _ in enumerate()
# print length
for megarep in range(MEGAREPS):
  for i in range(REPS):
    for temp in TEMPS:
      for time in TIMES:
        try: 
          # onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
          # onscale = int(onstring.split()[3])
          filename = "{0}{1}_{2}_{3}_{4}_on.tpkl".format(rootname, megarep+1, temp, i+1, time)
          on = parse(filename)
          # on_Nj = on.Nj
          on_SA, on_sigSA = on.scale_isosbestic()
          # print on.SA
          # print on_SA
          ons[time].append((on_SA, on_sigSA))
          # on_scaled = alg_scale(reference, on)
          filename = "{0}{1}_{2}_{3}_{4}_off.tpkl".format(rootname, megarep+1, temp, i+1, time)
          off = parse(filename)
          # on_Nj = on.Nj
          off_SA, off_sigSA = off.scale_isosbestic()
          offs[time].append((off_SA, off_sigSA))
        except:
          print "No data for Megarep {} and rep {}".format(megarep+1, i+1)

for time in TIMES:
  ### ONS
  ons_SA, ons_sigSA = zip(*ons[time])
  ons_SA_stdev = [np.std([i[j] for i in ons_SA]) for j in range(len(ons_SA[0]))]
  on_average_SA = sum(ons_SA)/len(ons_SA)
  # on_average_sigSA = np.power(sum([np.power(i, 2) for i in ons_sigSA])/len(ons_sigSA)**2,0.5)
  with open("sample_{}_on.dat".format(time), 'w') as file:
    for i, q in enumerate(on.q):
      file.write("{} {} {}\n".format(q, on_average_SA[i], ons_SA_stdev[i]))
  
  ### OFFS
  ons_SA, ons_sigSA = zip(*offs[time])
  ons_SA_stdev = [np.std([i[j] for i in ons_SA]) for j in range(len(ons_SA[0]))]
  on_average_SA = sum(ons_SA)/len(ons_SA)
  # on_average_sigSA = np.power(sum([np.power(i, 2) for i in ons_sigSA])/len(ons_sigSA),0.5)
  with open("sample_{}_off.dat".format(time), 'w') as file:
    for i, q in enumerate(on.q):
      file.write("{} {} {}\n".format(q, on_average_SA[i], ons_SA_stdev[i]))


    