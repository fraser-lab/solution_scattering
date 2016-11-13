#!/usr/bin/env python
'''
Simple buffer subtraction routine for PDH file generated
from Anton Paar SAXSess mc^2 machine.
The file skips the first 5 header lines of PDH file
and also throws out negative q values
The output is a "dat" file, a text file with three columns:
q, I(q), Error 
That can be directly used by GNOM for indirect fourier transformation

Usage:
subtractbuffer.py sample.pdh buffer.pdh output

04/04/2014 - Daniel (Agard Lab)
04/06/2014 - Added argument parser and buffer scaling switch
'''
import sys
import argparse
from math import sqrt

msg = '''
Simple buffer subtraction routine for PDH file generated
by Anton Paar SAXSess mc^2 machine.
The file skips the first 5 header lines of PDH file
and also throws out negative q values.
The output is a .dat file, a text file with three columns:
q,I(q),Error.
This can be directly used by GNOM for indirect fourier transformation\n
		'''
parser = argparse.ArgumentParser(description=msg, epilog="04/04/2014 - D.E.")
parser.add_argument("Sample", metavar="Sample.pdh", type=str, help="PDH file, sample scattering file.")
parser.add_argument("Buffer", metavar="Buffer.pdh", type=str, help="PDH file, buffer scattering file.")
parser.add_argument("Output", metavar="Output.dat", type=str, help="DAT file, buffer-subtracted sample")
parser.add_argument("-adj", "--adjustbuffer", action="store_true", help="A flag to adjust buffer scattering.")
parser.add_argument("-qmax", metavar="--qmax_cutoff", type=float, default=7.3, help="Cut off of high-angle data in nm^-1")
args = parser.parse_args()

def run(Sample_file, Buffer_file, Output_file, AdjustBuffer, Qmax):

	# Default buffer scaling is 1.0
	bufscale = 1.0

	s_fhd = open(Sample_file)
	b_fhd = open(Buffer_file)

	# read lines from input files
	s_raw = s_fhd.readlines()
	b_raw = b_fhd.readlines()

	# close file stream
	s_fhd.close()
	b_fhd.close()

	subtracted = [[], [], []]

	qlim_min = 0.05
	# For buffer scaling (if AdjustBuffer is true)
	qlim_bg_min = 6.5
	qlim_bg_max = 7.0

	# get raw q data for indexing
	rawq  = [float(s.strip('\r\n').split()[0]) for i,s in enumerate(s_raw) if i>4 and i<1285]

	qlim_min_id    = min([i for i,q in enumerate(rawq) if q>=qlim_min])
	qlim_bg_min_id = min([i for i,q in enumerate(rawq) if q>=qlim_bg_min])
	qlim_bg_max_id = min([i for i,q in enumerate(rawq) if q>=qlim_bg_max])
	qmax_id  	   = min([i for i,q in enumerate(rawq) if q>=Qmax])

	hiq_S = []
	hiq_B = []

	if AdjustBuffer:
		# calculate scaling factor
		# Assuming that q-spacing between sample and buffer is the same
		# We can just sum intensity values to approximate the integral
		# sum(S(q))/sum(B(q)), for 6.0 <= q <= 7.0 nm^-1
		for i in range(qlim_bg_min_id, qlim_bg_max_id+1):
			wrk_S = [float(d) for d in s_raw[i].strip('\r\n').split()]
			wrk_B = [float(d) for d in b_raw[i].strip('\r\n').split()]
			hiq_S.append(wrk_S[1])
			hiq_B.append(wrk_B[1])

		bufscale = sum(hiq_S)/sum(hiq_B)

		print "Calculated Buffer Scaling between 6 to 7 nm^-1.\nScale Factor : {0:.5f}".format(bufscale)

	q = []
	Iq = []
	Error = []

	for i in range(qlim_min_id,qmax_id):
		wrk_sample = [float(d) for d in s_raw[i].strip('\r\n').split()]
		wrk_buffer = [float(d) for d in b_raw[i].strip('\r\n').split()]
		q.append(wrk_sample[0])
		Iq.append(wrk_sample[1] - bufscale*wrk_buffer[1])
		Error.append(sqrt(wrk_sample[2]**2 + (bufscale*wrk_buffer[2])**2))

	out_fhd = open(Output_file,'wt')
	out_fhd.write('\n') # skip the first line

	for i in range(len(q)):
		rowout = "{0:1.9e} {1:1.9e} {2:1.9e}\n".format(q[i],Iq[i],Error[i])
		out_fhd.write(rowout)

	out_fhd.close()

	print "Subtraction Successful, wrote output to {0}.".format(Output_file)


if __name__=="__main__":
	Sample_file = args.Sample
	Buffer_file = args.Buffer
	Output_file = args.Output
	AdjustBuffer = args.adjustbuffer
	Qmax = args.qmax
	print "Qmax cut off is : {0:.2f}".format(Qmax)
	if '.dat' not in Output_file:
		Output_file = Output_file + '.dat'

	run(Sample_file, Buffer_file, Output_file, AdjustBuffer, Qmax)


