"""The 'trace' class is to work with traces from SAXS/WAXS data. 
It takes q and I or S and/or SA and provides simple utility classes 
for understanding trace data.

Benjamin Barad
"""



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
		print scale
		self.scaled_SA = self.SA/scale
		self.scaled_sigSA = self.sigSA/scale
		return (self.scaled_SA, self.scaled_sigSA)

	def as_vector(self):
		""" The SA column is the air-scattering adjusted integrated intensity"""
		return self.SA
	
	def get_q(self):
		return self.q
		
	def __repr__(self):
		final = "Q\tSA\tsigSA\tNj\n"
		for index, _ in enumerate(self.q):
			final += ("{0}\t{1}\t{2}\t{3}\n".format(self.q[index], self.SA[index],
						  self.sigSA[index],self.Nj[index]))
		return final
			
