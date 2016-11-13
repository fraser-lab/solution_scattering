import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import svd
import pandas
import csv




# sample_filename = "fast_times_sample_averaged.csv"
# sample_filename = "S99T_protein.csv"
sample_filename = "S99T_HD_outlier_free.csv"
sample_filename_2 = "WT_HD_protein_onoff_stdevs.csv"
# buffer_filename ="fast_times_buffer_averaged.csv"
# buffer_filename ="S99T_buffer.csv" #set_3_
buffer_filename = "S99T_HD_buffer_outlier_free.csv"
buffer_filename_2 = "WT_HD_buffer_outlier_free.csv"


# TIMES = ["10ns", "100ns", "1us", "10us"]

## 2/decade
# TIMES = ["-1us",  "316ns", "1us", "3.16us", "10us", "31.6us", "100us", "316us", "1ms",] # "10ns", "31.6ns", "100ns",  "3.16ms", "10ms"

## 4/decade
# TIMES = ["-1us",  "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us", "10us", "17.8us", "31.6us", "56.2us", "100us", "178us", "316us", "562us", "1ms", "1.78ms", "3.16ms", "5.62ms", "10ms"] # 
# TIMES = ["-10us",  "3.16us", "5.62us",  "10us", "17.8us", "31.6us", "56.2us", "100us", "316us", "562us"]
# TIMES = ["-1us", "10ns", "17.8ns", "31.6ns", "56.2ns", "75ns", "100ns", "133ns", "178ns", "316ns", "562ns", "1us", "1.78us", "3.16us", "5.62us" ]
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us",   "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us","422us", "562us", "750us", "1ms"] # 
# TIMES = ["-10.1us", "562ns", "1us", "1.78us",  "3.16us", "5.62us",  "10us",  "17.8us", "31.6us",  "56.2us",  "100us", "178us", "316us", "562us", "1ms"] # 

sample_table = pandas.read_csv(sample_filename, sep="\t", index_col="q")
sample_table_2 = pandas.read_csv(sample_filename_2, sep="\t", index_col="q")
buffer_table = pandas.read_csv(buffer_filename, sep="\t", index_col="q")
buffer_table_2 = pandas.read_csv(buffer_filename_2, sep="\t", index_col="q")
# print buffer_table[x]
q = sample_table.index.values
subtracted_dict = []
subtracted_dict_2 = []

fig, ax = plt.subplots()
plots=[]
for time in TIMES:
  print sample_table_2[time]
  subtracted_dict.append((time, (sample_table[time] - buffer_table[time]).values)) # 
  subtracted_dict_2.append((time, (sample_table_2[time] - buffer_table_2[time]).values))
  ax.plot(sample_table.index.values, [i*q for i,q in zip(((sample_table_2[time]-buffer_table_2[time])-(sample_table[time] - buffer_table[time])),sample_table.index.values)], label=time) # 
  #   ax.plot(sample_table.index.values, [i*q for i,q in zip(((sample_table_2[time]-buffer_table_2[time])-(sample_table[time] - buffer_table[time])),sample_table.index.values)], label=time) # 




ax.set_xscale("log", nonposx='clip')
ax.legend(loc='lower right', ncol=4, fontsize="small", framealpha=0.8)
ax.grid()
ax.set_ylabel(r"$\Delta\Delta\Delta I$")
ax.set_xlabel(r"Q ($nm^{-1}$)")
# plt.tight_layout()
# plt.show()
fig.savefig("WT-S99T-High-Outlers.png")


times_numeric = []
for time in TIMES[1:]:
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


fig2,ax2 = plt.subplots()
ax2.set_title("Integrated AUC", y=1.05)
ax2.plot(times_numeric, [-1*sum(subtracted_dict[i][1][0:12])*.0025 for i in range(1,len(TIMES))], ".", label="S99T AUC", color = "#5DA5DA")
ax2.plot(times_numeric, [-1*sum(subtracted_dict_2[i][1][0:12])*.0025 for i in range(1,len(TIMES))], ".", label="WT AUC", color = "#FAA43A")

# ax2.plot(times_numeric, [-1*sum(subtracted_dict[i][1][13:33]) for i in range(1,len(TIMES))], color="red", label="Integrated AUC from 0.05 to 0.1")
ax2.set_xlabel("Time (ns)")
ax2.set_xscale("log", nonposx='clip')
ax2.set_ylabel("AUC (q=0.02-0.04)")
ax2.grid()
ax2.legend(loc='upper left', ncol=1, framealpha=0.8)
# plt.tight_layout()
# plt.show()
fig2.savefig("Combined_HD_AUC_Plot_outliers.png")

fig3, ax3 = plt.subplots()
ax3.set_title("Difference in Integrated AUC between WT and S99T", y=1.05)
data_values = [-1*(sum(subtracted_dict_2[i][1][0:11])-sum(subtracted_dict[i][1][3:11]))*.0025 for i in range(1,len(TIMES))]
ax3.plot(times_numeric, data_values, ".", label="WT - S99T", color="#60BD68")
ax3.set_xlabel("Time (ns)")
ax3.set_xscale("log", nonposx='clip')
ax3.set_ylabel("Change in AUC(q=0.03-0.05)")
ax3.grid()
ax3.legend(loc='upper left', ncol=1, framealpha=0.8)
plt.tight_layout()
plt.show()
fig3.savefig("Difference_HD_AUC_Plot_outliers.png")

  
  
# with open("14C_total_subtraction.csv", 'wb') as csvfile:
#   writer = csv.writer(csvfile, delimiter="\t")
#   writer.writerow(["q"]+[i[0] for i in subtracted_dict])
#   for i in range(len(sample_table[time])):
#     # print i
#     writer.writerow([(i+86)*0.0025]+[j[1][i] for j in subtracted_dict])
#   