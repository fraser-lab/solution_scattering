from scipy.optimize import curve_fit
import numpy as np


def single_step_relaxation(x,a,b,c):
	# A is the asymptote
	# B is the kobs
	# C is the offset
	return a*(1-np.exp(-b*x))+c
	

def two_step_relaxation(x,a,b,c,d,e)
	# A and C are asymptotes
	# B and D are corresponding kobs
	# E is the floating offset
	return a*(1-np.exp(-b*x))+c*(1-np.exp(-d*x))+e

def three_step_relaxation(x,a,b,c,d,e,f,g)
	return a*(1-np.exp(-b*x))+c*(1-np.exp(-d*x))+e*(1-np.exp(-f*x))+g


def fit_kinetics(x,y):
	popt,pcov = curve_fit(single_step_relaxation,x,y, p0=[1,1,1], maxfev=5000)
	print popt
	print pcov
	kobs = popt[1]
	kobs_err = pcov[1,1]
	return popt,pcov