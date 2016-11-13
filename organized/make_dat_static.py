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
TEMPS = ["7", "14", "21", "28"]
# TIMES = ["-1us",  "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms", "1.78ms", "3.16ms", "5.62ms", "10ms"]
# TIMES = ["-1us"]
# MEGAREPS = 5
REPS = 28
MEGAREPS = 1
rootname = argv[1]

ons = {i: [] for i in TEMPS}

# files = listdir(directory)
# for index, _ in enumerate()
# print length
for megarep in range(MEGAREPS):
  for i in range(REPS):
    for temp in TEMPS:
        try: 
          # onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
          # onscale = int(onstring.split()[3])
          filename = "{0}{1}_{2}.tpkl".format(rootname, temp, i+1)
          on = parse(filename)
          # on_Nj = on.Nj
          on_SA, on_sigSA = on.scale_isosbestic()
          # print on.SA
          # print on_SA
          ons[temp].append((on_SA, on_sigSA))
          # on_scaled = alg_scale(reference, on)
        except:
          print "No data for Megarep {} and rep {}".format(megarep+1, i+1)

for temp in TEMPS:
  ### ONS
  ons_SA, ons_sigSA = zip(*ons[temp])
  ons_SA_stdev = [np.std([i[j] for i in ons_SA]) for j in range(len(ons_SA[0]))]
  # print ons_SA_stdev
  on_average_SA = sum(ons_SA)/len(ons_SA)
  on_average_sigSA = np.power(sum([np.power(i, 2) for i in ons_sigSA])/len(ons_sigSA)**2,0.5)
  print on_average_sigSA
  with open("static_buffer_{}C.dat".format(temp), 'w') as file:
    for i, q in enumerate(on.q):
      file.write("{} {} {}\n".format(q, on_average_SA[i], ons_SA_stdev[i]))
  
  # ### OFFS
  # ons_SA, ons_sigSA = zip(*offs[time])
  # ons_SA_stdev = [np.std([i[j] for i in ons_SA]) for j in range(len(ons_SA[0]))]
  # on_average_SA = sum(ons_SA)/len(ons_SA)
  # # on_average_sigSA = np.power(sum([np.power(i, 2) for i in ons_sigSA])/len(ons_sigSA),0.5)
  # with open("sample_{}_off.dat".format(time), 'w') as file:
  #   for i, q in enumerate(on.q):
  #     file.write("{} {} {}\n".format(q, on_average_SA[i], ons_SA_stdev[i]))


    