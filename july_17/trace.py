"""The 'trace' class is to work with traces from SAXS/WAXS data. 
It takes q and I or S and/or SA and provides simple utility classes 
for understanding trace data.

Benjamin Barad
"""

import numpy as np

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

	def alg_scale(self, ref, overwrite=None):
		"""Scale by projection of a vector onto a reference vector and determining the magnitude difference."""
		SA_ref = ref.SA
		SA_var = self.SA
		q = ref.q
		q_SA_ref = SA_ref*q
		q_SA_var = SA_var*q
		top = np.dot(q_SA_var,q_SA_ref)
		bottom = np.dot(q_SA_var,q_SA_var)
		scalar = top/bottom
		self.scale_factor = scalar
		if overwrite:
			self.SA = SA_var * scalar
			self.sigSA = self.sigSA * scalar
		else:
			self.scaled_SA = SA_var * scalar
			self.scaled_sigSA = self.sigSA * scalar
		return

	def integration_scale(self, ref):
		"""Scale by the total number of scattered photons"""
		SA_ref = ref.SA
		SA_var = self.SA
		q = self.q
		top = np.dot(SA_ref, q)
		bottom =  np.dot(SA_var, q)
		scalar = top/bottom
		self.scale_factor = scalar
		self.scaled_SA = SA_var * scalar
		self.scaled_sigSA = self.sigSA * scalar
		return

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
			
