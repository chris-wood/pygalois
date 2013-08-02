import pickle

class FF:
	''' Class for generic finite fields which can support arbitrary extensions
	and arithmetic in different bases.
	'''
	def __init__(self, base, exp = 1, ip = None):
		if (type(base) is int):

			# TODO: check if base is prime, otherwise not cyclic and can't be used to create the field

			self.isGroundField = True
			self.ground = base
		elif (type(base) is FF):
			# Test for irreducibility over base field
			# Irreducibility test: http://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Rabin.27s_test_of_irreducibility
			# Also check out HAC irreducibility test

			raise Exception("Not implemented")
		else:
			raise TypeError("Invalid base: " + str(base))

	def op_add(self, x, y):
		x = self.coerce(x)
		y = self.coerce(y)
		if (self.isGroundField):
			return (x + y) % self.ground
		else:
			raise Exception("TODO")

	def op_sub(self, x, y):
		x = self.coerce(x)
		y = self.coerce(y)
		if (self.isGroundField):
			return (x - y) % self.ground
		else:
			raise Exception("TODO")

	def op_mul(self, x, y):
		x = self.coerce(x)
		y = self.coerce(y)
		if (self.isGroundField):
			return (x * y) % self.ground
		else:
			raise Exception("TODO")

	def op_pow(self, x, k):
		x = self.coerce(x)
		if not (type(k) is int):
			raise Exception("Can only raise elements to integer powers.")
		if (self.isGroundField):
			return (x ** k) % self.ground
		else:
			raise Exception("TODO")

	def op_rem(self, x):
		x = self.coerce(x)
		if (self.isGroundField):
			return x % self.ground
		else:
			raise Exception("TODO")

	def op_div(self, x, y):
		x = self.coerce(x)
		y = self.coerce(y)
		if (self.isGroundField):
			yi = self.op_inv(y)
			return self.op_mul(x, yi)
		raise Exception("TODO")

	def op_inv(self, x):
		x = self.coerce(x)
		if (self.isGroundField):
			return (x ** (self.ground - 2)) % self.ground # Fermat's Little Theorem
		else:
			raise Exception("TODO")

	def coerce(self, x):
		if (self.isGroundField and type(x) is int):
			return x
		else:
			raise Exception("Cannot coerce element " + str(x) + " into field (ring) " + str(self))

	def eea(self, a, b):
		if (self.isGroundField):
			x = 0
			y = 1
			u = 1
			v = 0
			while b != 0:
				q = a / b
				a, b = b, a % b
				x, u = u - (q * x), x
				y, v = v - (q * y), y
			return u, v
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
		if (self.isGroundField):
			return None
		raise Exception("TODO")

	def getExtension(self):
		if (self.isGroundField):
			return 1
		raise Exception("TODO")

	def __str__(self):
		raise Exception("TODO")

