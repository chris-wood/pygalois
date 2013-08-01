import pickle

class FF:
	''' Class for generic finite fields which can support arbitrary extensions
	and arithmetic in different bases.
	'''
	def __init__(self, base, exp = 1, ip = None):
		if (type(base) is int):

			# TODO: implement isPrime() method

			self.isGroundField = True
			self.ground = base
		elif (type(base) is FF):
			# Test for irreducibility over base field

			raise TypeError("NOT DONE!")
			# Irreducibility test: http://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Rabin.27s_test_of_irreducibility

		raise Exception("TODO")

	def op_add(self, x, y):
		if (self.isGroundField):
			return (x + y) % self.ground
		raise Exception("TODO")

	def op_sub(self, x, y):
		if (self.isGroundField):
			return (x - y) % self.ground
		raise Exception("TODO")

	def op_mul(self, x, y):
		if (self.isGroundField):
			return (x * y) % self.ground
		raise Exception("TODO")

	def op_rem(self, x):
		if (self.isGroundField):
			return x % self.ground
		raise Exception("TODO")

	def op_div(self, x, y):
		raise Exception("TODO")

	def op_eea(self, a, b):
		raise Exception("TODO")

	def op_inv(self, x):
		raise Exception("TODO")

	def op_pow(self, x, k):
		if (self.isGroundField):
			return (x ** k) % self.ground
		raise Exception("TODO")

	def order(self, g): 
		raise Exception("TODO")

	def generators(self):
		raise Exception("TODO")

	def isGenerator(self, gen):
		raise Exception("TODO")

	def getFieldPolynomial(self):
		raise Exception("TODO")

	def getPrimitiveElement(self):
		raise Exception("TODO")

	def getSubfield(self):
		raise Exception("TODO")

	def getExtension(self):
		raise Exception("TODO")

	def __str__(self):
		raise Exception("TODO")