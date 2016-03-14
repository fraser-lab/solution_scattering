import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd
import pandas
import csv



sample_filename = "fast_times_sample_averaged.csv"
# sample_filename = "slow_times_sample_averaged.csv"
buffer_filename ="fast_times_buffer_averaged.csv"
# buffer_filename ="slow_times_buffer_averaged.csv"
TIMES = ["10ns", "100ns", "1us", "10us"]
# TIMES = ["-1us", "10us", "100us", "1ms", "10ms"]


sample_table = pandas.read_csv(sample_filename, sep="\t", index_col="x")
buffer_table = pandas.read_csv(buffer_filename, sep="\t", index_col="x")
# print buffer_table[x]
q = sample_table.index.values

subtracted_dict = []

fig, ax = plt.subplots()
plots=[]
for time in TIMES:
  subtracted_dict.append((time, (sample_table[time] - buffer_table[time]).values))
  plots.append(ax.plot(sample_table.index.values,sample_table[time] - buffer_table[time], label=time)[0])

ax.legend(plots)
plt.show()


  
with open("fast_times_total_subtraction.csv", 'wb') as csvfile:
  writer = csv.writer(csvfile, delimiter="\t")
  writer.writerow(["q"]+[i[0] for i in subtracted_dict])
  for i in range(len(sample_table[time])):
    print i
    writer.writerow([(i+6)*0.0025]+[j[1][i] for j in subtracted_dict])
  
