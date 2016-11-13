#! /usr/bin/env python


import math
import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl

# WT_filename = "WT_HD_protein.pkl"
WT_filename = "_data/WT_protein_chi.pkl"
with open(WT_filename, "rb") as pklfile:
  WT_protein = pkl.load(pklfile)

# Mut_filename = "S99T_HD_protein.pkl"
Mut_filename = "_data/S99T_proten_chi.pkl"
with open(Mut_filename, "rb") as pklfile:
  Mut_protein = pkl.load(pklfile)


# WT_buffer_filename = "WT_HD_buffer.pkl"
WT_buffer_filename = "_data/WT_buffer_chi.pkl"
with open(WT_buffer_filename, "rb") as pklfile:
  WT_buffer = pkl.load(pklfile)

# Mut_buffer_filename = "S99T_HD_buffer.pkl"
Mut_buffer_filename = "_data/S99T_buffer_chi.pkl"
with open(Mut_buffer_filename, "rb") as pklfile:
  Mut_buffer = pkl.load(pklfile)


q = [0.02+0.0025*i for i in range(2124)]
WT_subtracted = []
Mut_subtracted = []
double_subtracted = []

scale_factor = 0.956
mut_scale_factor = 0.92
variant = "WT"
time = "10us"
TIMES = ["-10.1us", "1us", "3.16us", "10us", "31.6us", "100us", "316us", "1ms"] # 
fig, ax = plt.subplots()

for index, value in enumerate(WT_protein):
  if value[0] == time:
    protein_means, protein_stds = zip(*value[1])
    buffer_means, buffer_stds = zip(*WT_buffer[index][1])

    ax.plot(q[4:], buffer_means[4:],  label = "Buffer: {}".format(time), color="#5DA5DA") #yerr=buffer_stds[4:],
    ax.plot(q[4:], protein_means[4:], label = "Protein: {}".format(time), color="#5DA5DA") #yerr=protein_stds[4:],
ax.set_title("Scattering Differences after T-Jump")
ax.legend(fontsize="x-small")
ax.set_xscale("log", nonposx='clip')
ax.set_ylabel(r"$\Delta$I (a.u.)")
ax.set_xlabel(r"q($\AA^{-1}$)")
ax.yaxis.set_ticks_position('none') # this one is optional but I still recommend it...
ax.xaxis.set_ticks_position('none')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_xlim(0.03,5.3275)
ax.set_xticks([0.03,0.1,1,5])
ax.set_xticklabels([0.03,0.1,1,5])
fig.savefig("On_off_difference.png")
plt.clf()
# plt.show()

# print WT_buffer[0][1]

fig, ax = plt.subplots()
ax.set_title("Buffer-subtracted Difference Scattering", y=1.05)
ax.set_xscale("log", nonposx='clip')
ax.legend(loc='lower right', ncol=4, fontsize="x-small", framealpha=0.8)
ax.set_ylabel(r"$\Delta\Delta I$ (a.u.)")
ax.set_xlabel(r"q ($\AA^{-1}$)")
ax.yaxis.set_ticks_position('none') # this one is optional but I still recommend it...
ax.xaxis.set_ticks_position('none')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_xlim(0.03,5.3275)
ax.set_ylim(-250, 50)
ax.set_xticks([0.03,0.1,1,5])
ax.set_xticklabels([0.03,0.1,1,5])

fig4, ax4 = plt.subplots()
fig5, ax5 = plt.subplots()

ax4.set_title("S99T Double Difference Scattering", y=1.05)
ax4.set_xscale("log", nonposx='clip')
ax4.set_ylabel(r"$\Delta\Delta I$ (a.u.)")
ax4.set_xlabel(r"q ($\AA^{-1}$)")
ax4.yaxis.set_ticks_position('none') # this one is optional but I still recommend it...
ax4.xaxis.set_ticks_position('none')
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.get_xaxis().tick_bottom()
ax4.get_yaxis().tick_left()
ax4.set_xlim(0.03,5.3275)
ax4.set_ylim(-250, 50)
ax4.set_xticks([0.03,0.1,1,5])
ax4.set_xticklabels([0.03,0.1,1,5])

ax5.set_title("Difference between S99T and WT", y=1.05)
ax5.set_xscale("log", nonposx='clip')
ax5.legend(loc='lower right', ncol=4, fontsize="x-small", framealpha=0.8)
ax5.set_ylabel(r"$\Delta\Delta\Delta I$ (a.u.)")
ax5.set_xlabel(r"q ($\AA^{-1}$)")
ax5.yaxis.set_ticks_position('none') # this one is optional but I still recommend it...
ax5.xaxis.set_ticks_position('none')
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.get_xaxis().tick_bottom()
ax5.get_yaxis().tick_left()
ax5.set_xlim(0.03,5.3275)
ax5.set_xticks([0.03,0.1,1,5])
ax5.set_xticklabels([0.03,0.1,1,5])

fig6, ax6 = plt.subplots()
ax6.set_title("Buffer-subtracted Difference Scattering", y=1.05)
ax6.set_xscale("log", nonposx='clip')
ax6.legend(loc='lower right', ncol=4, fontsize="x-small", framealpha=0.8)
ax6.set_ylabel(r"$\Delta\Delta I$ (a.u.)")
ax6.set_xlabel(r"q ($\AA^{-1}$)")
ax6.yaxis.set_ticks_position('none') # this one is optional but I still recommend it...
ax6.xaxis.set_ticks_position('none')
ax6.spines['top'].set_visible(False)
ax6.spines['right'].set_visible(False)
ax6.get_xaxis().tick_bottom()
ax6.get_yaxis().tick_left()
ax6.set_xlim(0.03,5.3275)
ax6.set_ylim(-250, 50)
ax6.set_xticks([0.03,0.1,1,5])
ax6.set_xticklabels([0.03,0.1,1,5])


plots=[]
for index, value in enumerate(WT_protein):
  WT_means = [WT_protein[index][1][k][0]-(WT_buffer[index][1][k][0]*scale_factor) for k in range(len(WT_protein[index][1]))]
  WT_stds = [np.sqrt(WT_protein[index][1][k][1]**2+(WT_buffer[index][1][k][1]*scale_factor)**2) for k in range(len(WT_protein[index][1]))]
  WT_subtracted.append((value[0], WT_means, WT_stds)) # 
  ##
  Mut_means = [Mut_protein[index][1][k][0]-(Mut_buffer[index][1][k][0]*mut_scale_factor) for k in range(len(Mut_protein[index][1]))]
  Mut_stds = [np.sqrt(Mut_protein[index][1][k][1]**2+(Mut_buffer[index][1][k][1]*mut_scale_factor)**2) for k in range(len(Mut_protein[index][1]))]
  Mut_subtracted.append((value[0], Mut_means, Mut_stds)) 
  ## 
  subtracted_means = [WT_means[i] - Mut_means[i] for i,_ in enumerate(Mut_means)]
  subtracted_stds = [np.sqrt(Mut_stds[i]**2+WT_stds[i]**2) for i,_ in enumerate(Mut_stds)]
  double_subtracted.append((value[0], subtracted_means, subtracted_stds))
  if value[0] == time:
    ax6.plot(q, WT_means, label=value[0], color = "#5DA5DA")
    ax6.legend(fontsize="x-small")
    fig6.savefig("Double_subtraction_one_time.png")
  if value[0] in TIMES:
    ax.plot(q, WT_means,  label=value[0]) #yerr=WT_stds,
    ax4.plot(q, Mut_means,  label=value[0]) #yerr=WT_stds,
    ax5.plot(q, subtracted_means,  label=value[0]) #yerr=WT_stds,
  # ax4.plot(q, Mut_means,  label=value[0]) #yerr=Mut_stds,
  # ax5.plot(q, subtracted_means,  label=value[0]) #yerr=subtracted_stds,
  # ax.plot(q, WT_means, label=value[0])
  # ax.plot(sample_table.index.values, [i*q for i,q in zip(((sample_table_2[time]-buffer_table_2[time])-(sample_table[time] - buffer_table[time])),sample_table.index.values)], label=time) # 
ax.legend(loc='lower right', ncol=4, fontsize="x-small", framealpha=0.8)
ax4.legend(loc='lower right', ncol=4, fontsize="x-small", framealpha=0.8)
ax5.legend(loc='lower right', ncol=4, fontsize="x-small", framealpha=0.8)
fig.savefig("WT_double_differences.png")
plt.clf()
fig4.savefig("Mut_double_differences.png")
fig5.savefig("Triple_differences.png")





# ax4.set_title("Mutant")
# ax4.set_xscale("log", nonposx='clip')
# ax4.legend(loc='lower right', ncol=4, fontsize="small", framealpha=0.8)
# ax4.grid()
# ax4.set_ylabel(r"$\Delta\Delta\Delta I$")
# ax4.set_xlabel(r"Q ($nm^{-1}$)")

# ax5.set_title("Subtracted")
# ax5.set_xscale("log", nonposx='clip')
# ax5.legend(loc='lower right', ncol=4, fontsize="small", framealpha=0.8)
# ax5.grid()
# ax5.set_ylabel(r"$\Delta\Delta\Delta I$")
# ax5.set_xlabel(r"Q ($nm^{-1}$)")
# # plt.tight_layout()
# # plt.show()
# fig.savefig("3.png")


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

print times_numeric

fig2,ax2 = plt.subplots()
ax2.set_title("Integrated Area of Guinier Region", y=1.05)
ax2.set_xlabel("Time (ns)")
ax2.set_xscale("log", nonposx='clip')
ax2.set_ylabel("Area (q=0.03-0.06)")
ax2.yaxis.set_ticks_position('none') # this one is optional but I still recommend it...
ax2.xaxis.set_ticks_position('none')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
ax2.set_ylim(2.5, 5.5)
ax2.set_xlim(500, 2000000)

fig7,ax7 = plt.subplots()
ax7.set_title("Integrated Area of Guinier Region", y=1.05)
ax7.set_xlabel("Time (ns)")
ax7.set_xscale("log", nonposx='clip')
ax7.set_ylabel("Area (q=0.03-0.06)")
ax7.yaxis.set_ticks_position('none') # this one is optional but I still recommend it...
ax7.xaxis.set_ticks_position('none')
ax7.spines['top'].set_visible(False)
ax7.spines['right'].set_visible(False)
ax7.get_xaxis().tick_bottom()
ax7.get_yaxis().tick_left()
ax7.set_ylim(2.5,5.5)
ax7.set_xlim(500,2000000)


# integrated_AUCs = [-1*sum(WT)]
WT_integrated = [-1*sum(WT_subtracted[i][1][4:17])*.0025 for i in range(1,len(times_numeric))]
WT_integrated_errors=[-1* np.sqrt(sum([k**2 for k in WT_subtracted[i][2][4:17]]))*.0025 for i in range(1,len(times_numeric))]
WT = ax2.errorbar(times_numeric[1:], WT_integrated, fmt=".", yerr=WT_integrated_errors, label="Integrated Scattering", color = "#5DA5DA")
ax2.legend(loc='upper left', fontsize="x-small", ncol=1, framealpha=0.8)
fig2.savefig("WT_integrated.png")
print "---"
print WT_integrated
print "---"
Mut_integrated = [-1*sum(Mut_subtracted[i][1][4:17])*.0025 for i in range(1,len(times_numeric))]
Mut_integrated_errors=[-1* np.sqrt(sum([k**2 for k in Mut_subtracted[i][2][4:17]]))*.0025 for i in range(1,len(times_numeric))]
print Mut_integrated
ax7.errorbar(times_numeric[1:], Mut_integrated, fmt=".", yerr=Mut_integrated_errors, label="S99T Integrated Scattering", color = "#FAA43A")
ax7.legend(loc='upper left', fontsize="x-small", ncol=1, framealpha=0.8)
fig7.savefig("S99T_integrated.png")

ax7.errorbar(times_numeric[1:], WT_integrated, fmt=".", yerr=Mut_integrated_errors, label="WT Integrated Scattering", color = "#5DA5DA")
ax7.legend(loc='upper left', fontsize="x-small", ncol=1, framealpha=0.8)

# # ax2.plot(times_numeric, [-1*sum(subtracted_dict[i][1][13:33]) for i in range(1,len(TIMES))], color="red", label="Integrated AUC from 0.05 to 0.1")

# # plt.tight_layout()
# # plt.show()
fig7.savefig("Paired diference.png")

# x = range(1000,1000000)

fig3, ax3 = plt.subplots()
ax3.set_title("Integrated WT-S99T Guinier Difference", y=1.05)
double_integrated = [-1*sum(double_subtracted[i][1][4:17])*.0025 for i in range(1,len(times_numeric))]
double_integrated_errors=[-1* np.sqrt(sum([k**2 for k in double_subtracted[i][2][2:17]]))*.0025 for i in range(1,len(times_numeric))]
ax3.errorbar(times_numeric[1:], double_integrated, fmt=".", yerr=double_integrated_errors, label="WT-S99T Integrated Scattering", color = "#5DA5DA")
# ax3.plot(x, np.poly1d(np.polyfit(times_numeric[, y, 1))(x))
ax3.set_xlabel("Time (ns)")
ax3.set_xscale("log", nonposx='clip')
ax3.set_ylabel("Change in Area (q=0.03-0.06)")
ax3.legend(loc='upper left', fontsize="x-small", ncol=1, framealpha=0.8)
# ax3.set_xlim(500,2000000)
ax3.yaxis.set_ticks_position('none') # this one is optional but I still recommend it...
ax3.xaxis.set_ticks_position('none')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.get_xaxis().tick_bottom()
ax3.get_yaxis().tick_left()
# plt.tight_layout()
# # plt.show()
fig3.savefig("difference_integrated.png")


DI_max = np.mean(WT_integrated[-6:-1])
DI_t = [max(WT_integrated)-i for i in WT_integrated]
DI_0 = max(WT_integrated) - min(WT_integrated)
y = [i/DI_0 for i in DI_t]
x = times_numeric[1:]
fig, ax = plt.subplots()
x = x[11:19]
y = [math.log(i) for i in y[11:19]]
m,b = np.polyfit(x,y,1)
print m, b
ax.plot(x, y, ".")
ax.plot(x, [m*i + b for i in x], "-")
ax.set_title("Exponential Decay Fits")
# ax.set_yscale("log", nonposy="clip")
ax.set_ylabel("ln(Fraction Remaining/Total Decay)")
ax.set_xlabel("Time (ns)")
ax.yaxis.set_ticks_position('none') # this one is optional but I still recommend it...
ax.xaxis.set_ticks_position('none')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
fig.savefig("WT_best_fit.png")
plt.show()
  
# with open("14C_total_subtraction.csv", 'wb') as csvfile:
#   writer = csv.writer(csvfile, delimiter="\t")
#   writer.writerow(["q"]+[i[0] for i in subtracted_dict])
#   for i in range(len(sample_table[time])):
#     # print i
#     writer.writerow([(i+86)*0.0025]+[j[1][i] for j in subtracted_dict])
#   