"""The 'trace' class is to work with traces from SAXS/WAXS data. 
It takes q and I or S and/or SA and provides simple utility classes 
for understanding trace data.

Benjamin Barad
"""

import numpy as np
from numpy import linalg as LA

class Trace(object):
	def __init__(self, q, sigS, S, sigSA, SA, Nj):
		self.q = q
		self.sigS = sigS
		self.S = S
		self.sigSA = sigSA
		self. SA = SA
		self.Nj = Nj

	def scale(self, ref, qmin=0.0025, qmax=5.1925, approach="algebraic"):
		"""
                Scale by the total number of scattered photons
		
		:param ref: Reference curve (other Trace object)
        :param qmin: Minimum q for scaling (default 0.0025) 
		:param qmax: Maximum q for scaling (default 5.1925)
        :returns: Scaled SA and sigSA curves (which also get written on the trace as self.scaled_SA, self.scaled_sigSA)
		"""
		SA_var = self.SA[self.q>=qmin] ###0.0025 works for full q across all cypa
		SA_var = SA_var[self.q<=qmax]  ###5.1925 works for full q across all cypa
		SA_ref = ref.SA[self.q>=qmin]
		SA_ref = ref.SA[self.q<=qmax]

		if approach == "algebraic":
			q = ref.q[self.q>=qmin]
			q = ref.q[self.q<=qmax]
			q_SA_ref = SA_ref*q
			q_SA_var = SA_var*q
			top = np.dot(q_SA_var,q_SA_ref)
			bottom = np.dot(q_SA_var,q_SA_var)
			scalar = top/bottom

		elif approach == "projection":
			ref_norm = LA.norm(SA_ref)
			ref_hat = SA_ref / ref_norm
			scalar_projection = np.dot(SA_var, ref_hat)
			scalar = scalar_projection / ref_norm

		elif approach == "integration":
			Nj = ref.Nj[self.q>=qmin]
			Nj = ref.Nj[self.q<=qmax]
			top = np.dot(SA_var,Nj)
			bottom = np.dot(SA_ref,Nj)
			scalar = top/bottom

		else:
			raise TypeError('please choose from the following scaling approaches: algebraic, projection, integration')

		self.scale_factor = scalar
		self.scaled_SA = self.SA * scalar
		self.scaled_sigSA = self.sigSA * scalar
		return self.scaled_SA, self.scaled_sigSA

	def buffer_scale(self, ref):
		"""Scale by the total number of scattered photons"""
		SA_ref = ref.SA
		SA_var = self.SA
		q = self.q
		data_mask = np.array(q, dtype=bool)
		data_mask[q<1.5]=False
		data_mask[q>3.6]=False
		q = q[data_mask]
		SA_ref = SA_ref[data_mask]
		SA_var = SA_var[data_mask]
		top = np.dot(SA_ref, q)
		bottom =  np.dot(SA_var, q)
		scalar = top/bottom
		self.buffer_scale_factor = scalar
		self.scaled_SA = self.SA * scalar
		self.scaled_sigSA = self.sigSA * scalar
		return

	def subtract(self, trace_two, scaled=None, buffer_scaled=None):
		if scaled:
			err_one = (self.scale_factor*self.sigSA)**2
			err_two = (trace_two.scale_factor*trace_two.sigSA)**2
			err_cov = (2*self.scale_factor*trace_two.scale_factor*np.cov(self.sigSA,trace_two.sigSA)[0][1])
			output_SA = (self.scaled_SA - trace_two.scaled_SA)
		elif buffer_scaled:
			err_one = (self.buffer_scale_factor*self.sigSA)**2
			err_two = (trace_two.buffer_scale_factor*trace_two.sigSA)**2
			err_cov = (2*self.buffer_scale_factor*trace_two.buffer_scale_factor*np.cov(self.sigSA,trace_two.sigSA)[0][1])
			output_SA = (self.scaled_SA - trace_two.scaled_SA)
		else:
			err_one = self.sigSA**2
			err_two = trace_two.sigSA**2
			err_cov = (2*np.cov(self.sigSA,trace_two.sigSA)[0][1])
			output_SA = (self.SA - trace_two.SA)
		total_err = np.sqrt(np.abs(err_one + err_two - err_cov))
		output = Trace(self.q, np.empty_like(self.q), np.empty_like(self.q), total_err, output_SA, np.empty_like(self.q))

		return output


	def add(self, trace_two, scaled=None, buffer_scaled=None):
		if scaled:
			err_one = (self.scale_factor*self.sigSA)**2
			err_two = (trace_two.scale_factor*trace_two.sigSA)**2
			err_cov = (2*self.scale_factor*trace_two.scale_factor*np.cov(self.sigSA,trace_two.sigSA)[0][1])
			output_SA = (self.scaled_SA + trace_two.scaled_SA)
		elif buffer_scaled:
			err_one = (self.buffer_scale_factor*self.sigSA)**2
			err_two = (trace_two.buffer_scale_factor*trace_two.sigSA)**2
			err_cov = (2*self.buffer_scale_factor*trace_two.buffer_scale_factor*np.cov(self.sigSA,trace_two.sigSA)[0][1])
			output_SA = (self.scaled_SA + trace_two.scaled_SA)
		else:
			err_one = self.sigSA**2
			err_two = trace_two.sigSA**2
			err_cov = (2*np.cov(self.sigSA,trace_two.sigSA)[0][1])
			output_SA = (self.SA + trace_two.SA)
		total_err = np.sqrt(np.abs(err_one + err_two + err_cov))
		output = Trace(self.q, np.empty_like(self.q), np.empty_like(self.q), total_err, output_SA, np.empty_like(self.q))

		return output

	def as_vector(self):
		""" The SA column is the air-scattering adjusted integrated intensity"""
		return self.SA
	
	def get_q(self):
		return self.q

	def set_name(self, name):
		self.name=name
		return

	def write_dat(self, filename):
	    data = np.column_stack((self.q, self.SA, self.sigSA))
	    np.savetxt(filename, data, fmt='%f', delimiter='    ', newline='\n', header='q    I    sigI')
	    print("data successfully written to {}".format(filename))
	    return
		
	def __repr__(self):
		final = "Q\tSA\tsigSA\tNj\n"
		for index, _ in enumerate(self.q):
			final += ("{0}\t{1}\t{2}\t{3}\n".format(self.q[index], self.SA[index],
						  self.sigSA[index],self.Nj[index]))
		return final

