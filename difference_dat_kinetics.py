# Generate kinetics from time-resolved difference dat files
import csv

import numpy as np
from matplotlib import pyplot as plt

from relax import relaxation_fit, fit_bootstrap_two, fit_bootstrap_single, single_step_relaxation, two_step_relaxation, three_step_relaxation
from parse import parse
import trace

FOLDER = "CypA-WT-1_March"
PREFIX = "CypA-WT-1_diff_" #
# PREFIX = "_data/"
TIMES_STR = ["-10.1us", "562ns","750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us","23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us"] # "750us", "1ms"
#INITIAL_GUESS = (-1, 1./1000, 1, 1./10000, 1)
#INITIAL_GUESS = (-1, 1./1000, 1)
INITIAL_GUESS = (.8,1./10000, 2)
 # 1./1000000, 3) # (-1, 1./10000, 1)(.8,1./10000, 2) #
RELAXATION_STEPCOUNT = single_step_relaxation #two_step_relaxation #single_step_relaxation 
QMIN = 0.03
QMAX = 0.05


def time_str_to_float(time_string):
  number = float(time_string[:-2])
  scale = time_string[-2:]
  if scale == "ns":
    scaled_number = number
  elif scale == "us":
    scaled_number = 1000 * number
  elif scale == "ms":
    scaled_number = 1000 * 1000 * number
  else:
    print("scale could not be calculated")
  return scaled_number


def plot_integrated_areas(tuple_list, y_calc, y_calc_2, filename = "integrated_area_over_time.png"):
	# Tuple list should be of the form [(time_numeric, trace, integrated area, integrated error),(),...]
	fig, ax = plt.subplots()
	x, _, y, yerr = zip(*tuple_list)
	curve = ax.errorbar(x[1:],y[1:], fmt=".", yerr=yerr[1:])
	curve_2 = ax.plot(x[1:], y_calc, '-')
	# curve_3 = ax.plot(x[1:], y_calc_2, '-')
	ax.set_xscale('log')
	ax.set_xlim(x[1], x[-1])
	fig.savefig(filename)
	return fig,ax

def plot_differences(tuple_list, filename="differences.png"):
	fig, ax = plt.subplots()
	_, traces, _, _ = zip(*tuple_list)
	for index, trace in enumerate(traces):
		ax.plot(trace.q, trace.SA, "-", label=TIMES_STR[index])
	ax.legend()
	ax.set_xscale('log')
	ax.set_xlim(np.min(traces[0].q), np.max(traces[0].q))
	fig.savefig(filename)
	return fig,ax


def integrate_area(trace, q_min = QMIN, q_max = QMAX):
	q = trace.get_q()
	index_low = np.nonzero(q>=q_min)[0][0]
	index_high = np.nonzero(q<=q_max)[0][-1]

	series_I = []
	series_error = []
	for i in range(index_low, index_high+1):
		delta_q = q[i+1] - q[i]
		# print(delta_q)
		I_section = trace.SA[i] * delta_q
		error_section = trace.sigSA[i] * delta_q
		series_I.append(I_section)
		series_error.append(error_section)
	integrated_area = -1*sum(series_I)
	integrated_error = np.sqrt(sum([i**2 for i in series_error]))
	return integrated_area, integrated_error


def measure_kinetics(area_series, time_series, initial, funct=RELAXATION_STEPCOUNT, maxfev=30000, sigma=None):
	x = time_series
	y = area_series
	popt, pcov, y_calc = relaxation_fit(x,y, relaxation_function=funct, initial_guess=initial, maxfev=maxfev, sigma=sigma)
	return popt, pcov, y_calc

def write_csv(times, areas, errors):
	assert len(times)==len(areas)==len(errors)
	with open("integrated_areas_modified_wt_single_new.csv", "w") as csvfile:
		writer = csv.writer(csvfile, delimiter=",")
		for i,time in enumerate(times):
			writer.writerow([time, areas[i], errors[i]])


def run(prefix, times_str):
	traces = []

	for time in times_str:
		trace = parse("{2}/{0}{1}.dat".format(prefix, time, FOLDER))
		time_numeric = time_str_to_float(time)
		area, error = integrate_area(trace)
		traces.append((time_numeric,trace, area, error))
	times,_,areas, errors = zip(*traces)
	write_csv(times, areas, errors)
	parameters, covariances, y_calc = measure_kinetics(areas[1:], times[1:], initial=INITIAL_GUESS, funct=RELAXATION_STEPCOUNT, sigma=errors[1:])
	print("Parameters of Fit:")
	if RELAXATION_STEPCOUNT == two_step_relaxation:
		print("First Step:")
		print("A1: {}\tStandard Deviation: {}".format(parameters[0], np.sqrt(covariances[0][0])))
		print("kobs1: {}\tStandard Deviation: {}".format(parameters[1], np.sqrt(covariances[1][1])))
		print("A2: {}\tStandard Deviation: {}".format(parameters[2], np.sqrt(covariances[2][2]))) 
		print("kobs2: {}\tStandard Deviation: {}".format(parameters[3], np.sqrt(covariances[3][3])))
		print("offset: {}\tStandard Deviation: {}".format(parameters[4], np.sqrt(covariances[4][4])))
		bootstrap_fit, bootstrap_error, y_calc_2 = fit_bootstrap_two(areas[1:], times[1:], initial_guess=parameters, sigma=errors[1:])

		# print (bootstrap_fit)
		# print (bootstrap_error)
		print("Bootstrap:")
		print("A1: {}\tStandard Deviation: {}".format(bootstrap_fit[0], bootstrap_error[0]))
		print("kobs1: {}\tStandard Deviation: {}".format(bootstrap_fit[1], bootstrap_error[1]))
		print("A2: {}\tStandard Deviation: {}".format(bootstrap_fit[2], bootstrap_error[2])) 
		print("kobs2: {}\tStandard Deviation: {}".format(bootstrap_fit[3], bootstrap_error[3]))
		print("offset: {}\tStandard Deviation: {}".format(bootstrap_fit[4], bootstrap_error[4]))

	if RELAXATION_STEPCOUNT == single_step_relaxation:
		print("First Step:")
		print("A1: {}\tStandard Deviation: {}".format(parameters[0], np.sqrt(covariances[0][0])))
		print("kobs1: {}\tStandard Deviation: {}".format(parameters[1], np.sqrt(covariances[1][1])))
		print("offset: {}\tStandard Deviation: {}".format(parameters[2], np.sqrt(covariances[2][2])))
		bootstrap_fit, bootstrap_error, y_calc_2 = fit_bootstrap_single(areas[1:], times[1:], initial_guess=parameters, sigma=errors[1:])

		# print (bootstrap_fit)
		# print (bootstrap_error)
		print("Bootstrap:")
		print("A1: {}\tStandard Deviation: {}".format(bootstrap_fit[0], bootstrap_error[0]))
		print("kobs1: {}\tStandard Deviation: {}".format(bootstrap_fit[1], bootstrap_error[1]))
		print("offset: {}\tStandard Deviation: {}".format(bootstrap_fit[2], bootstrap_error[2]))

	else:
		print(parameters)





	#print("Errors Sum/Areas Sum: {}".format(sum(errors)/sum(areas)))
	#print("K2 STD/ K2obs: {}".format(np.sqrt(covariances[3][3])/parameters[3]))
	plot_integrated_areas(traces, y_calc, y_calc_2)
	plot_differences(traces)


if __name__ == "__main__":
	run(PREFIX, TIMES_STR)

