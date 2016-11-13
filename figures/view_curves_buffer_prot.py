from sys import argv

import pandas as pd

import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
import numpy as np

variant = "WT"
time = "10us"

protein_off = pd.read_csv("_data/datfiles/{}_protein_off.dat".format(variant), sep=" ", header=None, names=["q", "I", "sigI"])
buffer_off = pd.read_csv("_data/datfiles/{}_buffer_off.dat".format(variant), sep=" ", header=None, names=["q", "I", "sigI"])

protein_on = pd.read_csv("_data/datfiles/{}_protein_{}_on.dat".format(variant, time), sep=" ", header=None, names=["q", "I", "sigI"])
buffer_on = pd.read_csv("_data/datfiles/{}_buffer_{}_on.dat".format(variant, time), sep=" ", header=None, names=["q", "I", "sigI"])


fig, ax = plt.subplots()
fig.set_size_inches(4,3)
ax.set_xscale("log", nonposx='clip')
buffer_off_plot, = ax.plot(buffer_off.q[4:], buffer_off.I[4:], label="Buffer: Dark")
buffer_on_plot, = ax.plot(buffer_on.q[4:], buffer_on.I[4:], label = "Buffer: {}".format(time))

ax.set_xlim(1.2,3.5)
ax.set_ylim(1500,4000)
ax.yaxis.set_ticks_position('none') # this one is optional but I still recommend it...
ax.xaxis.set_ticks_position('none')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xticks([1.5,2,2.5,3])
ax.set_xticklabels([1.5,2,2.5,3])
ax.set_yticks([2000,3000,4000])
ax.set_yticklabels([2000,3000,4000])
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

fig.savefig("Buffer_curves_inset.png")
fig.set_size_inches(8,6)
xlabel = ax.set_ylabel("Intensity (a.u.)")
ylabel = ax.set_xlabel(r"q ($\AA^{-1}$)")
ax.set_xticks([0.03,0.1,1,5])
ax.set_xticklabels([0.03,0.1,1,5])
ax.set_yticks(range(0,10001,2000))
ax.set_yticklabels(range(0,10001,2000))


ax.legend(fontsize='x-small')
ax.yaxis.set_ticks_position('left') # this one is optional but I still recommend it...
ax.xaxis.set_ticks_position('bottom')
ax.set_xlim(0.03, 5.3)
ax.set_ylim(0, 10000)
ax.legend(fontsize='x-small')
fig.savefig("Buffer_curves.png")
rectangle = ax.add_patch(matplotlib.patches.Rectangle((1.2, 1500), 2.3, 2500, fill=None))
fig.savefig("Buffer_curves_rectangle.png")
rectangle.remove()


# plt.show()
protein_off_plot, = ax.plot(protein_off.q[4:], protein_off.I[4:], label="CypA: Dark")
protein_on_plot, = ax.plot(protein_on.q[4:], protein_on.I[4:], label="CypA: {}".format(time))

ax.legend(fontsize='x-small')
fig.savefig("Protein_Buffer_curves.png")
rectangle = ax.add_patch(matplotlib.patches.Rectangle((0.03, 7000), 0.03, 3000, fill=None))
fig.savefig("Protein_curves_rectangle.png")
rectangle.remove()

ax.set_xticks([0.03,0.04,.05])
ax.set_xticklabels([0.03,0.04,0.05])
ax.set_yticks(range(7000,10001,1000))
ax.set_yticklabels(range(7000,10001,1000))

ax.set_xlabel("")
ax.set_ylabel("")
ax.legend_.remove()
ax.set_xlim(0.03, 0.06)
ax.set_ylim(7000,10000)
fig.set_size_inches(4,3)
fig.savefig("Protein_curves_inset.png")
# plt.show()