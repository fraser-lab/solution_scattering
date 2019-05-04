from scipy.optimize import curve_fit
import string
from math import sqrt
from matplotlib import pyplot as plt
from numpy import linspace, log

recip_T_1_list = [0.003350084, 0.003461405, 0.003516174, 0.003369272, 0.00333667, 0.003448276, 0.003395586, 0.00330142]
k1_list = [533273.0, 489773.0, 431018.0, 1165650.0, 977200.0, 583069.0, 572787.0, 1157360.0]
k1_err_list = [329011.0, 96143.0, 74712.0, 331986.0, 116888.0, 108298.0, 70106.0, 240739.0]

ln_of_k1_over_T_list = []
ln_of_k1_over_T_err_list = []

for i, rate1 in enumerate(k1_list):
  ln_of_k1_over_T = log(rate1*recip_T_1_list[i])
  ln_of_k1_over_T_list.append(ln_of_k1_over_T)
  ln_of_k1_over_T_err = k1_err_list[i] / k1_list[i]
  ln_of_k1_over_T_err_list.append(ln_of_k1_over_T_err)
  

recip_T_2_list = [0.003350084, 0.003461405, 0.003516174, 0.003369272, 0.00333667, 0.003448276, 0.003395586, 0.00330142, 0.003581662]
k2_list = [26555.0, 24936.0, 11856.0, 22317.0, 25403.0, 15603.0, 21728.0, 23779.0, 10964.0]
k2_err_list = [8678.0, 2586.0, 1761.0, 2251.0, 1636.0, 1616.0, 1729.0, 2383.0, 2857.0]

ln_of_k2_over_T_list = []
ln_of_k2_over_T_err_list = []

for j, rate2 in enumerate(k2_list):
  ln_of_k2_over_T = log(rate2*recip_T_2_list[i])
  ln_of_k2_over_T_list.append(ln_of_k2_over_T)
  ln_of_k2_over_T_err = k2_err_list[i] / k2_list[i]
  ln_of_k2_over_T_err_list.append(ln_of_k2_over_T_err)


#filename = "k1_Eyring.png"

filename = "k1_k2_Eyring_05032019_bootstrap_final.png"

def line( x, a, b ):
  return (a*x)+b

params1, covar1 = curve_fit(line, recip_T_1_list, ln_of_k1_over_T_list, sigma=ln_of_k1_over_T_err_list, absolute_sigma=True)
params2, covar2 = curve_fit(line, recip_T_2_list, ln_of_k2_over_T_list, sigma=ln_of_k2_over_T_err_list, absolute_sigma=True)

#params, covar = curve_fit(line, recip_T, ln_of_k_over_T)

DH_over_R_1 = 0.0-float(params1[0])
SD_DH_over_R_1 = sqrt(float(covar1[0,0]))

DH_1 = DH_over_R_1*8.3144598
SD_DH_1 = SD_DH_over_R_1*8.3144598

DS_1 = 8.3144598*(float(params1[1])-19.15480729345)
SD_DS_1 = 8.3144598*(sqrt(float(covar1[1,1])))

print "The enthalpy of activation for k1 is: %s +/- %s" % (DH_1, SD_DH_1)
print "The entropy of activation for k1 is: %s +/- %s" % (DS_1, SD_DS_1)

DH_over_R_2 = 0.0-float(params2[0])
SD_DH_over_R_2 = sqrt(float(covar2[0,0]))

DH_2 = DH_over_R_2*8.3144598
SD_DH_2 = SD_DH_over_R_2*8.3144598

DS_2 = 8.3144598*(float(params2[1])-19.15480729345)
SD_DS_2 = 8.3144598*(sqrt(float(covar2[1,1])))

print "The enthalpy of activation for k2 is: %s +/- %s" % (DH_2, SD_DH_2)
print "The entropy of activation for k2 is: %s +/- %s" % (DS_2, SD_DS_2)

x1 = []
y1 = ln_of_k1_over_T_list
err1 = ln_of_k1_over_T_err_list
for value1 in recip_T_1_list:
  value1 = 1000.0*value1
  x1.append(value1)

x2 = []
y2 = ln_of_k2_over_T_list
err2 = ln_of_k2_over_T_err_list
for value2 in recip_T_2_list:
  value2 = 1000.0*value2
  x2.append(value2)

fig, ax = plt.subplots(figsize=(6,8))
ax.scatter(x1,y1, s=350, marker='.', color='black', zorder=3)
linspace_min1 = 1000*recip_T_1_list[-1]
linspace_max1 = 1000*recip_T_1_list[2]
fit_x1_vals = linspace(linspace_min1, linspace_max1, num=100)
ax.plot(fit_x1_vals, [(params1[0]/1000.0)*i+params1[1] for i in fit_x1_vals], "-", color='#9970ab', linewidth='3',zorder=1)
ax.scatter(x2,y2, s=350, marker='.', color='black', zorder=3)
linspace_min2 = 1000*recip_T_2_list[-2]
linspace_max2 = 1000*recip_T_2_list[-1]
fit_x2_vals = linspace(linspace_min2, linspace_max2, num=100)
ax.plot(fit_x2_vals, [(params2[0]/1000.0)*j+params2[1] for j in fit_x2_vals], "-", color='#5aae61', linewidth='3',zorder=1)
ax.set_xlabel("1000*1/T")
ax.set_ylabel("ln(k/T)")
ax.errorbar(x1, y1, fmt=".", yerr=err1, color='black', linewidth='2', zorder=2)
ax.errorbar(x2, y2, fmt=".", yerr=err2, color='black', linewidth='2', zorder=2)
ax.yaxis.set_tick_params(width=2)
ax.yaxis.set_ticks_position('left') # this one is optional but I still recommend it...
ax.xaxis.set_tick_params(width=2)
ax.xaxis.set_ticks_position('bottom')
plt.xlim(3.28, 3.62)
plt.ylim(3.0, 9.0)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
  item.set_fontsize(16)#ax.set_ylabel("Ylabel Goes Here", labelpad=10)
#ax.set_xlabel("XLabel goes here",labelpad=10)
plt.tight_layout()
fig.savefig(filename)

