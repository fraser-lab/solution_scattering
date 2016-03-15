import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd
import pandas
import csv




# sample_filename = "fast_times_sample_averaged.csv"
sample_filename = "set_6_protein_averaged_scaled_isosbestic.csv"
# buffer_filename ="fast_times_buffer_averaged.csv"
buffer_filename ="set_6_buffer_averaged_scaled_isosbestic.csv" #set_3_


# TIMES = ["10ns", "100ns", "1us", "10us"]

## 2/decade
# TIMES = ["-1us",  "316ns", "1us", "3.16us", "10us", "31.6us", "100us", "316us", "1ms",] # "10ns", "31.6ns", "100ns",  "3.16ms", "10ms"

## 4/decade
# TIMES = ["-1us",  "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms", "1.78ms", "3.16ms", "5.62ms", "10ms"] # 
# TIMES = ["-1us", "10us", "17.8us", "31.6us", "56.2us", "100us", "562us"]
TIMES =["-1us", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms", "3.16ms"]

times_numeric = []
for time in TIMES:
	number = float(time[:-2])
	scale = time[-2:]
	if scale == "ns":
		times_numeric.append(number)
	elif scale == "us":
		times_numeric.append(1000*number)
	elif scale == "ms":
		times_numeric.append(1000*1000*number)
	else:
		print "scale could not be calculated"


# TIMES = ["-1us", "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us" ]
sample_table = pandas.read_csv(sample_filename, sep="\t", index_col="q")
buffer_table = pandas.read_csv(buffer_filename, sep="\t", index_col="q")
# print buffer_table[x]
q = sample_table.index.values
subtracted_dict = []

fig, ax = plt.subplots()
plots=[]
for time in TIMES:
  scale = sum([abs(sample_table[time][i]*i) for i in sample_table.index.values if i>0.3])/sum([abs(buffer_table[time][i]*i) for i in sample_table.index.values if i>0.3])
  print time, scale
  subtracted_dict.append((time, (sample_table[time] - buffer_table[time]*scale).values)) # 
  ax.plot(sample_table.index.values, sample_table[time] - buffer_table[time]*scale , label=time) # 

ax.set_xscale("log", nonposx='clip')
ax.set_title("Buffer-subtracted on-off differences", y=1.05)
ax.legend(loc='lower right', ncol=2, framealpha=0.8)
ax.grid()
ax.set_ylabel(r"$\Delta\Delta I$")
ax.set_xlabel(r"Q ($\AA ^{-1}$)")
ax.yaxis.tick_left()
ax.xaxis.tick_bottom()
plt.tight_layout()
# plt.show(block=False)
fig.savefig("14C_all.png")

fig2,ax2 = plt.subplots()
ax2.set_title("Integrated AUC", y=1.05)
ax2.plot(times_numeric, [-1*sum(subtracted_dict[i][1][1:13]) for i in range(len(TIMES))], label="Integrated AUC from 0.02 to 0.05")
ax2.plot(times_numeric, [-1*sum(subtracted_dict[i][1][13:33]) for i in range(len(TIMES))], color="red", label="Integrated AUC from 0.05 to 0.1")
ax2.set_xlabel("Time (ns)")
ax2.set_xscale("log", nonposx='clip')
ax2.set_ylabel("Negative AUC from 0.02 to 0.1")
ax2.grid()
# ax2.set_ylim(4,6)
ax2.legend(loc='lower left', ncol=1, framealpha=0.8)
plt.tight_layout()
# plt.show()
fig2.savefig("AUC_Plot.png")



  
# with open("14C_total_subtraction.csv", 'wb') as csvfile:
#   writer = csv.writer(csvfile, delimiter="\t")
#   writer.writerow(["q"]+[i[0] for i in subtracted_dict])
#   for i in range(len(sample_table[time])):
#     # print i
#     writer.writerow([(i+86)*0.0025]+[j[1][i] for j in subtracted_dict])
#   