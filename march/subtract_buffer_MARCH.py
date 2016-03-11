import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd
import pandas
import csv




# sample_filename = "fast_times_sample_averaged.csv"
sample_filename = "set_6_protein_averaged_scaled.csv"
# buffer_filename ="fast_times_buffer_averaged.csv"
buffer_filename ="set_6_buffer_averaged_scaled.csv" #set_3_


# TIMES = ["10ns", "100ns", "1us", "10us"]

## 2/decade
# TIMES = ["-1us",  "316ns", "1us", "3.16us", "10us", "31.6us", "100us", "316us", "1ms",] # "10ns", "31.6ns", "100ns",  "3.16ms", "10ms"

## 4/decade
# TIMES = ["-1us",  "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms", "1.78ms", "3.16ms", "5.62ms", "10ms"] # 
TIMES = ["-1us", "10us", "17.8us", "31.6us", "56.2us", "100us", "562us"]
# TIMES = ["-1us", "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us" ]
sample_table = pandas.read_csv(sample_filename, sep="\t", index_col="q")
buffer_table = pandas.read_csv(buffer_filename, sep="\t", index_col="q")
# print buffer_table[x]
q = sample_table.index.values
subtracted_dict = []

fig, ax = plt.subplots()
plots=[]
for time in TIMES:
  subtracted_dict.append((time, (sample_table[time] - buffer_table[time] ).values)) # 
  ax.plot(sample_table.index.values, sample_table[time] - buffer_table[time] , label=time) # 

ax.set_xscale("log", nonposx='clip')
ax.legend(loc='lower right', ncol=2, framealpha=0.8)
ax.grid()
ax.set_ylabel(r"$\Delta I$")
ax.set_xlabel(r"Q ($nm^{-1}$)")
plt.tight_layout()
plt.show()
fig.savefig("14C_all.png")


  
# with open("14C_total_subtraction.csv", 'wb') as csvfile:
#   writer = csv.writer(csvfile, delimiter="\t")
#   writer.writerow(["q"]+[i[0] for i in subtracted_dict])
#   for i in range(len(sample_table[time])):
#     # print i
#     writer.writerow([(i+86)*0.0025]+[j[1][i] for j in subtracted_dict])
#   