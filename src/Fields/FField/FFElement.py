class FFElement(object):
	''' A class to store a generic element in GF(p^n), where p is prime or some
	smaller field. Elements in GF(p) for p prime are represented as integers and
	not as objects of this type, so there is no need to do any extra heap allocation.
	'''
	def __init__(self, field, coeffs):
		self.coeffs = coeffs[:]
		self.field = field
