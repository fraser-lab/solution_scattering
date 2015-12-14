"""The 'trace' class is to work with traces from SAXS/WAXS data. 
It takes q and I or S and/or SA and provides simple utility classes 
for understanding trace data.

Benjamin Barad
"""

class Trace(object):
	def __init__(self, q, sigS, S, sigSA, SA):
		self.q = q
		self.sigS = sigS
		self.S = S
		self.sigSA = sigSA
		self. SA = SA

	def apply_i0_scaling(self, scale_factor):
		self.scaled_SA = self.SA/scale_factor
		return self.scaled_SA

	def as_vector(self):
		""" The SA column is the air-scattering adjusted integrated intensity"""
		return self.SA
	
	def get_q(self):
		return self.q
		
	def __repr__(self):
		final = ""
		for index, _ in enumerate(self.q):
			final += ("{0}\t{1}\t{2}\n".format(self.q[index], self.S[index],
						  self.SA[index]))
		return final
			
