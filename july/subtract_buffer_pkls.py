import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl

# WT_filename = "WT_HD_protein.pkl"
WT_filename = "WT_protein.pkl"
with open(WT_filename, "rb") as pklfile:
  WT_protein = pkl.load(pklfile)

# Mut_filename = "S99T_HD_protein.pkl"
Mut_filename = "S99T_protein.pkl"
with open(Mut_filename, "rb") as pklfile:
  Mut_protein = pkl.load(pklfile)


# WT_buffer_filename = "WT_HD_buffer.pkl"
WT_buffer_filename = "WT_buffer.pkl"
with open(WT_buffer_filename, "rb") as pklfile:
  WT_buffer = pkl.load(pklfile)


# Mut_buffer_filename = "S99T_HD_buffer.pkl"
Mut_buffer_filename = "S99T_buffer.pkl"
with open(Mut_buffer_filename, "rb") as pklfile:
  Mut_buffer = pkl.load(pklfile)


q = [0.02+0.0025*i for i in range(2124)]
WT_subtracted = []
Mut_subtracted = []
double_subtracted = []

# print WT_buffer[0][1]

fig, ax = plt.subplots()
fig4, ax4 = plt.subplots()
fig5, ax5 = plt.subplots()
plots=[]
for index, value in enumerate(WT_protein):
  WT_means = [WT_protein[index][1][k][0]-WT_buffer[index][1][k][0] for k in range(len(WT_protein[index][1]))]
  WT_stds = [np.sqrt(WT_protein[index][1][k][1]**2+WT_buffer[index][1][k][1]**2) for k in range(len(WT_protein[index][1]))]
  WT_subtracted.append((value[0], WT_means, WT_stds)) # 
  ##
  Mut_means = [Mut_protein[index][1][k][0]-Mut_buffer[index][1][k][0] for k in range(len(Mut_protein[index][1]))]
  Mut_stds = [np.sqrt(Mut_protein[index][1][k][1]**2+Mut_buffer[index][1][k][1]**2) for k in range(len(Mut_protein[index][1]))]
  Mut_subtracted.append((value[0], Mut_means, Mut_stds)) 
  ## 
  subtracted_means = [WT_means[i] - Mut_means[i] for i,_ in enumerate(Mut_means)]
  subtracted_stds = [np.sqrt(Mut_stds[i]**2+WT_stds[i]**2) for i,_ in enumerate(Mut_stds)]
  double_subtracted.append((value[0], subtracted_means, subtracted_stds))
  
  ax.errorbar(q, WT_means, yerr=WT_stds, label=value[0]) #
  ax4.errorbar(q, Mut_means, yerr=Mut_stds, label=value[0]) 
  ax5.errorbar(q, subtracted_means, yerr=subtracted_stds, label=value[0])
  # ax.plot(q, WT_means, label=value[0])
  # ax.plot(sample_table.index.values, [i*q for i,q in zip(((sample_table_2[time]-buffer_table_2[time])-(sample_table[time] - buffer_table[time])),sample_table.index.values)], label=time) # 



ax.set_title("WT")
ax.set_xscale("log", nonposx='clip')
ax.legend(loc='lower right', ncol=4, fontsize="small", framealpha=0.8)
ax.grid()
ax.set_ylabel(r"$\Delta\Delta\Delta I$")
ax.set_xlabel(r"Q ($nm^{-1}$)")

ax4.set_title("Mutant")
ax4.set_xscale("log", nonposx='clip')
ax4.legend(loc='lower right', ncol=4, fontsize="small", framealpha=0.8)
ax4.grid()
ax4.set_ylabel(r"$\Delta\Delta\Delta I$")
ax4.set_xlabel(r"Q ($nm^{-1}$)")

ax5.set_title("Subtracted")
ax5.set_xscale("log", nonposx='clip')
ax5.legend(loc='lower right', ncol=4, fontsize="small", framealpha=0.8)
ax5.grid()
ax5.set_ylabel(r"$\Delta\Delta\Delta I$")
ax5.set_xlabel(r"Q ($nm^{-1}$)")
# plt.tight_layout()
# plt.show()
fig.savefig("3.png")


times_numeric = []
for time in zip(*WT_protein)[0]:
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
# integrated_AUCs = [-1*sum(WT)]
WT_integrated = [-1*sum(WT_subtracted[i][1][4:17])*.0025 for i in range(1,len(times_numeric))]
WT_integrated_errors=[-1* np.sqrt(sum([k**2 for k in WT_subtracted[i][2][4:17]]))*.0025 for i in range(1,len(times_numeric))]
ax2.errorbar(times_numeric[1:], WT_integrated, fmt=".", yerr=WT_integrated_errors, label="WT AUC", color = "#5DA5DA")

Mut_integrated = [-1*sum(Mut_subtracted[i][1][4:17])*.0025 for i in range(1,len(times_numeric))]
Mut_integrated_errors=[-1* np.sqrt(sum([k**2 for k in Mut_subtracted[i][2][4:17]]))*.0025 for i in range(1,len(times_numeric))]
ax2.errorbar(times_numeric[1:], Mut_integrated, fmt=".", yerr=Mut_integrated_errors, label="S99T AUC", color = "#FAA43A")

# # ax2.plot(times_numeric, [-1*sum(subtracted_dict[i][1][13:33]) for i in range(1,len(TIMES))], color="red", label="Integrated AUC from 0.05 to 0.1")
ax2.set_xlabel("Time (ns)")
ax2.set_xscale("log", nonposx='clip')
ax2.set_ylabel("AUC (q=0.025-0.06)")
ax2.grid()
ax2.legend(loc='upper left', ncol=1, framealpha=0.8)
# # plt.tight_layout()
# # plt.show()
fig2.savefig("1.png")

x = range(1000,1000000)

fig3, ax3 = plt.subplots()
ax3.set_title("Difference in Integrated AUC between WT and S99T", y=1.05)
double_integrated = [-1*sum(double_subtracted[i][1][4:17])*.0025 for i in range(1,len(times_numeric))]
double_integrated_errors=[-1* np.sqrt(sum([k**2 for k in double_subtracted[i][2][4:17]]))*.0025 for i in range(1,len(times_numeric))]
ax3.errorbar(times_numeric[1:], double_integrated, fmt=".", yerr=double_integrated_errors, label="Difference AUC", color = "#FAA43A")
# ax3.plot(x, np.poly1d(np.polyfit(times_numeric[, y, 1))(x))
ax3.set_xlabel("Time (ns)")
ax3.set_xscale("log", nonposx='clip')
ax3.set_ylabel("Change in AUC(q=0.03-0.06)")
ax3.grid()
ax3.legend(loc='upper left', ncol=1, framealpha=0.8)

plt.tight_layout()
# plt.show()
fig3.savefig("2.png")
  
  
# with open("14C_total_subtraction.csv", 'wb') as csvfile:
#   writer = csv.writer(csvfile, delimiter="\t")
#   writer.writerow(["q"]+[i[0] for i in subtracted_dict])
#   for i in range(len(sample_table[time])):
#     # print i
#     writer.writerow([(i+86)*0.0025]+[j[1][i] for j in subtracted_dict])
#   