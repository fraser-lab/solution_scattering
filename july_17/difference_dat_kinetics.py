# Generate kinetics from time-resolved difference dat files

import numpy as np
from matplotlib import pyplot as plt

from parse import parse
import trace

PREFIX = "_data/CypA-WT-1_diff_"
TIMES_STR = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us","23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"]


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


def plot_integrated_areas(tuple_list, filename = "integrated_area_over_time.png"):
	# Tuple list should be of the form [(time_numeric, trace, integrated area, integrated error),(),...]
	fig, ax = plt.subplots()
	x, _, y, yerr = zip(*tuple_list)
	curve = ax.errorbar(x[1:],[-i for i in y[1:]], fmt=".", yerr=yerr[1:])
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


def integrate_area(trace, q_min = 0.03, q_max = 0.06):
	q = trace.get_q()
	index_low = np.nonzero(q == q_min)[0][0]
	index_high = np.nonzero(q==q_max)[0][0]

	series_I = []
	series_error = []
	for i in range(index_low, index_high+1):
		delta_q = q[i+1] - q[i]
		I_section = trace.SA[i] * delta_q
		error_section = trace.sigSA[i] * delta_q
		series_I.append(I_section)
		series_error.append(error_section)
	integrated_area = sum(series_I)
	integrated_error = np.sqrt(sum([i**2 for i in series_error]))
	return integrated_area, integrated_error


def measure_kinetics(trace):
	return trace

def run(prefix, times_str):
	traces = []
	for time in times_str:
		trace = parse("{0}{1}.dat".format(prefix, time))
		time_numeric = time_str_to_float(time)
		area, error = integrate_area(trace)
		traces.append((time_numeric,trace, area, error))
	plot_integrated_areas(traces)
	plot_differences(traces)


if __name__ == "__main__":
	run(PREFIX, TIMES_STR)

