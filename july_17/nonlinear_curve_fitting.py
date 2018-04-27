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


def 