from scipy.optimize import curve_fit
import string
from math import sqrt
from matplotlib import pyplot as plt

recip_T = [0.003350084, 0.003461405, 0.003581662, 0.003516174, 0.003369272, 0.00333667, 0.003448276, 0.003395586, 0.00330142]
ln_of_k_over_T = [4.4103, 4.3382, 3.4947, 3.7191, 4.3026, 4.4840, 3.9489, 4.3178, 4.4161]
ln_of_k_over_T_err = [0.165391191, 0.148223487930567, 0.291399369359574, 0.150106609808102, 0.0980803428936209, 0.0792813829987571, 0.108348843392715, 0.0929622086444897, 0.0846328720137199]

filename = "k2_Eyring.png"

#filename = "k2_Eyring_noSigmaInFit.png"

def line( x, a, b ):
  return (a*x)+b

params, covar = curve_fit(line, recip_T, ln_of_k_over_T, sigma=ln_of_k_over_T_err)

#params, covar = curve_fit(line, recip_T, ln_of_k_over_T)


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
