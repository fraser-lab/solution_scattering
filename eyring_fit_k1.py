from scipy.optimize import curve_fit
import string
from math import sqrt
from matplotlib import pyplot as plt

recip_T = [0.003350084, 0.003461405, 0.003516174, 0.003369272, 0.00333667, 0.003448276, 0.003395586, 0.00330142]
ln_of_k_over_T = [7.794164785, 7.358590555, 7.23359504, 8.179644904, 8.016032644, 7.602706005, 7.521864765, 8.282137606]
ln_of_k_over_T_err = [0.743780384, 0.351342113, 0.363095223, 0.343034962, 0.214385481, 0.264613706, 0.184722153, 0.193365145]

#filename = "k1_Eyring.png"

filename = "k1_Eyring_noSigmaInFit.png"

def line( x, a, b ):
  return (a*x)+b

#params, covar = curve_fit(line, recip_T, ln_of_k_over_T, sigma=ln_of_k_over_T_err)

params, covar = curve_fit(line, recip_T, ln_of_k_over_T)

DH_over_R = 0.0-float(params[0])
SD_DH_over_R = sqrt(float(covar[0,0]))

DH = DH_over_R*8.3144598
SD_DH = SD_DH_over_R*8.3144598

DS = 8.3144598*(float(params[1])-19.15480729345)
SD_DS = 8.3144598*(sqrt(float(covar[1,1])))

print "The enthalpy of activation is: %s +/- %s" % (DH, SD_DH)
print "The entropy of activation is: %s +/- %s" % (DS, SD_DS)

x = []
y = ln_of_k_over_T
err = ln_of_k_over_T_err
for value in recip_T:
  value = 1000.0*value
  x.append(value)

fig, ax = plt.subplots()
ax.plot(x,y, ".")
ax.plot(x, [(params[0]/1000.0)*i+params[1] for i in x], ls="--")
ax.set_xlabel("1000*1/T")
ax.set_ylabel("ln(k/T)")
ax.errorbar(x, y, fmt=".", yerr=err)
fig.savefig(filename)

