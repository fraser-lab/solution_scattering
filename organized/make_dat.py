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
TEMPS = ["14C", "14.001C"]
TIMES = ["-1us",  "10ns", "31.6ns", "100ns", "316ns", "1us", "3.16us", "10us", "31.6us", "100us", "316us", "1ms", "3.16ms", "10ms"] # 
# TIMES = ["-1us"]
MEGAREPS = 5
REPS = 5
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
  ons_SA, ons_sigSA = zip(*ons[time])
  on_average_SA = sum(ons_SA)/len(ons_SA)
  on_average_sigSA = np.power(sum([np.power(i, 2) for i in ons_sigSA])/len(ons_sigSA)**2,0.5)
  with open("buffer_{}_on.dat".format(time), 'w') as file:
    for i, q in enumerate(on.q):
      file.write("{} {} {}\n".format(q, on_average_SA[i], on_average_sigSA[i]))
  ons_SA, ons_sigSA = zip(*offs[time])
  on_average_SA = sum(ons_SA)/len(ons_SA)
  on_average_sigSA = np.power(sum([np.power(i, 2) for i in ons_sigSA])/len(ons_sigSA)**2,0.5)
  with open("buffer_{}_off.dat".format(time), 'w') as file:
    for i, q in enumerate(on.q):
      file.write("{} {} {}\n".format(q, on_average_SA[i], on_average_sigSA[i]))
    