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

	def apply_i0_scaling(self, scale_factor):
		self.scaled_SA = self.SA/scale_factor
		return self.scaled_SA

	def scale_isosbestic(self):
		Q = [0.0175 + 0.0025 * i for i in range(2125)]
		scale = sum([self.SA[i]*Q[i] for i in range(583,603)])
		print(scale)
		self.scaled_SA = self.SA/scale
		self.scaled_sigSA = self.sigSA/scale
		return (self.scaled_SA, self.scaled_sigSA)

	def alg_scale(self, ref, qmin=0.0025, qmax=5.1925):
		"""
                Scale by projection of a vector onto a reference vector and determining the magnitude difference.
		
		:param ref: Reference curve (other Trace object)
                :param qmin: Minimum q for scaling (default 0.0025) 
		:param qmax: Maximum q for scaling (default 5.1925)
                :returns: Scaled SA and sigSA curves (which also get written on the trace as self.scaled_SA, self.scaled_sigSA)
		"""
		SA_var = self.SA[self.q>=qmin] ###0.0025 works for full q across all cypa
		SA_var = SA_var[self.q<=qmax]  ###5.1925 works for full q across all cypa
		q = ref.q[self.q>=qmin]
		q = ref.q[self.q<=qmax]
		SA_ref = ref.SA[self.q>=qmin]
		SA_ref = ref.SA[self.q<=qmax]
		q_SA_ref = SA_ref*q
		q_SA_var = SA_var*q
		top = np.dot(q_SA_var,q_SA_ref)
		bottom = np.dot(q_SA_var,q_SA_var)
		scalar = top/bottom
		self.scale_factor = scalar
		self.scaled_SA = self.SA * scalar
		self.scaled_sigSA = self.sigSA * scalar
		return self.scaled_SA, self.scaled_sigSA

	def projection_scale(self, ref, qmin=0.0025, qmax=5.1925):
		"""
                Scale factor is determined by the ratio of - the scalar projection of a vector onto a reference vector - over - the norm of the reference vector
		
		:param ref: Reference curve (other Trace object)
                :param qmin: Minimum q for scaling (default 0.0025) 
		:param qmax: Maximum q for scaling (default 5.1925)
                :returns: Scaled SA and sigSA curves (which also get written on the trace as self.scaled_SA, self.scaled_sigSA)
		"""

		SA_var = self.SA[self.q>=qmin] ###0.0025 works for full q across all cypa
		SA_var = SA_var[self.q<=qmax]  ###5.1925 works for full q across all cypa
		SA_ref = ref.SA[self.q>=qmin]
		SA_ref = ref.SA[self.q<=qmax]

		ref_norm = LA.norm(SA_ref)
		ref_hat = SA_ref / ref_norm
		scalar_projection = np.dot(SA_var, ref_hat)
		scalar = scalar_projection / ref_norm
		self.scale_factor = scalar
		self.scaled_SA = self.SA * scalar
		self.scaled_sigSA = self.sigSA * scalar
		return

	def integration_scale(self, ref, qmin=0.0025, qmax=5.1925):
		"""
                Scale by the total number of scattered photons
		
		:param ref: Reference curve (other Trace object)
                :param qmin: Minimum q for scaling (default 0.0025) 
		:param qmax: Maximum q for scaling (default 5.1925)
                :returns: Scaled SA and sigSA curves (which also get written on the trace as self.scaled_SA, self.scaled_sigSA)
		"""
		SA_var = self.SA[self.q>=qmin] ###0.0025 works for full q across all cypa
		SA_var = SA_var[self.q<=qmax]  ###5.1925 works for full q across all cypa
		Nj = ref.Nj[self.q>=qmin]
		Nj = ref.Nj[self.q<=qmax]
		SA_ref = ref.SA[self.q>=qmin]
		SA_ref = ref.SA[self.q<=qmax]
		top = np.dot(SA_var,Nj)
		bottom = np.dot(SA_ref,Nj)
		scalar = top/bottom
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

	def subtract(self, trace_two):
	    err_one = self.sigSA**2
	    err_two = trace_two.sigSA**2
	    err_cov = (2*np.cov(self.sigSA,trace_two.sigSA)[0][1])
	    total_err = np.sqrt(np.abs(err_one+err_two-err_cov))
	    output_SA = (self.SA - trace_two.SA)
	    output = Trace(self.q, np.empty_like(self.q), np.empty_like(self.q), total_err, output_SA, np.empty_like(self.q))

	    err_one = (trace_one.scale_factor*trace_one.sigSA)**2
	    err_two = (trace_two.scale_factor*trace_two.sigSA)**2
	    err_cov = (2*trace_one.scale_factor*trace_two.scale_factor*np.cov(trace_one.sigSA,trace_two.sigSA)[0][1])
	    total_err = np.sqrt(np.abs(err_one+err_two-err_cov))
	    output_SA = (trace_one.scaled_SA - trace_two.scaled_SA)
	    output = Trace(trace_one.q, np.empty_like(trace_one.q), np.empty_like(trace_one.q), total_err, output_SA, np.empty_like(trace_one.q))
	    return output

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

class TraceMethods(Trace):
	def __init__(self, q=None, sigS=None, S=None, sigSA=None, SA=None, Nj=None):
		super().__init__(q, sigS, S, sigSA, SA, Nj)
		self.q = q
		self.sigS = sigS
		self.S = S
		self.sigSA = sigSA
		self. SA = SA
		self.Nj = Nj
