"""The 'trace' class is to work with traces from SAXS/WAXS data. 
It takes q and I or S and/or SA and provides simple utility classes 
for understanding trace data.

Benjamin Barad
"""

class Trace(object):
	def __init__(self, q, S, SA):
		self.q = q
		self.S = S
		self. SA = SA
		
	def as_vector(self):
		return self.SA
		
	def __repr__(self):
		final = ""
		for index, _ in enumerate(self.q):
			final += ("{0}\t{1}\t{2}\n".format(self.q[index], self.S[index],
						  self.SA[index]))
		return final
			
