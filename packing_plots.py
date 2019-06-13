from matplotlib import pyplot as plt

packing_dat_1 = 'CypA-WT-static-1_8C_spf.dat'
packing_dat_2 = 'CypA-WT-static-1_13C_spf.dat'
packing_dat_3 = 'CypA-WT-static-1_18C_spf.dat'
packing_dat_4 = 'CypA-WT-static-1_23C_spf.dat'
packing_dat_5 = 'CypA-WT-static-1_28C_spf.dat'

packing_dat_list =[packing_dat_1, packing_dat_2, packing_dat_3, packing_dat_4, packing_dat_5]

q_list = []
packing_list_1 = []
packing_list_2 = []
packing_list_3 = []
packing_list_4 = []
packing_list_5 = []

packing_lists = [packing_list_1, packing_list_2, packing_list_3, packing_list_4, packing_list_5]

temp_list = [8.0, 13.0, 18.0, 23.0, 28.0]
A2_list = []

for i, dat_file in enumerate(packing_dat_list):
  file = open(dat_file, 'r')
  for line in file.readlines(3:):
    if i==0:
      q_val = line.split()[0]
      q_list.append(q_val)
      intensity = line.split()[1]
      packing_list_1.append
    else:
      intensity = line.split()[1]
      packing_lists[i].append(intensity)

fig, ax = plt.subplots(figsize=(6,6))
ax.plot(q_list, packing_list_1, "-", color='#feb24c', linewidth='2',zorder=0)
ax.plot(q_list, packing_list_2, "-", color='#fd8d3c', linewidth='2',zorder=0)
ax.plot(q_list, packing_list_3, "-", color='#fc4e2a', linewidth='2',zorder=0)
ax.plot(q_list, packing_list_4, "-", color='#e31a1c', linewidth='2',zorder=0)
ax.plot(q_list, packing_list_5, "-", color='#bd0026', linewidth='2',zorder=0)

ax.yaxis.set_tick_params(width=2)
ax.yaxis.set_ticks_position('left') # this one is optional but I still recommend it...
ax.xaxis.set_tick_params(width=2)
ax.xaxis.set_ticks_position('bottom')
ax.set_xscale('log')
#plt.xlim(0.03, 10)
#plt.ylim(0, 10000)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
  ax.spines[axis].set_zorder(5)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
  item.set_fontsize(16)#ax.set_ylabel("Ylabel Goes Here", labelpad=10)
#ax.set_xlabel("XLabel goes here",labelpad=10)
plt.tight_layout()
fig.savefig("packing_structure_factor_wrt_T.png")





