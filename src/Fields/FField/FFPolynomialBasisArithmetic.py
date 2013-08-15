import pickle
import sys

class FFPolynomialBasisArithmetic:
	''' Class that implements polynomial basis arithmetic for finite fields.
	'''
	def __init__(self, field):
		raise Exception("TODO")

	def g_add(self, x, y):
		sum = []
		index = 0
		for i in range(len(x) + len(y)):
			z = (x[i] + y[i]) % self.base
			sum.insert(0, z)
			index = index + 1
		result = FFElement(sum)
		return result

	def g_sub(self, x, y):
		diff = []
		index = 0
		for i in range(len(x) + len(y)):
			z = (x[i] - y[i]) % self.base
			diff.insert(0, z)
			index = index + 1
		return FFElement(diff)

	# def g_mult(self, x, y, reduce = True):
	# 	prod = []
	# 	for i in range(len(x) + len(y)):
	# 		prod.append(0)
	# 	for i in range(len(x)):
	# 		for j in range(len(y)): 
	# 			pl = i + j
	# 			tempProd = 0
	# 			if (i < len(x) and j < len(y)):
	# 				tempProd = x[i] * y[j]
	# 			pl = i + j
	# 			prod[pl] = prod[pl] + tempProd
	# 	prod.reverse()

	# 	# Reduce the elements not of the highest degree
	# 	foundNonZero = False
	# 	if (reduce == True):
	# 		for i in range(len(prod)):
	# 			prod[i] = prod[i] % self.base

	# 	result = FFElement(prod)
	# 	if (result.degree() >= self.exp and reduce == True):
	# 		result = self.g_div(result, self.ip)[1] # div returns (Q,R), we want R(emainder)
	# 	return result

	# def reduce(self, x):
	# 	result = x.copy()
	# 	if (result.degree() >= self.exp):
	# 		result = self.g_div(result, self.ip)[1] # div returns (Q,R), we want R(emainder)
	# 	return result

	# def g_div(self, N, D):
	# 	if (D.degree() < 0):
	# 		raise TypeError()
	# 	if (N.degree() >= D.degree()): # numerator degree larger than divisor
	# 		Q = GFElem([]) # q <- 0
	# 		while (N.degree() >= D.degree()):
	# 			divisor = D.shiftLeft(N.degree() - D.degree()) # find the divisor for the degrees - this is correct
	# 			Q[N.degree() - D.degree()] = N[N.degree()] / divisor[divisor.degree()]
	# 			divisor = divisor.scalarMult(Q[N.degree() - D.degree()]) 
	# 			N = self.g_sub(N, divisor)
	# 		R = N # left over numerator is the remainder...
	# 		return (Q, R)
	# 	else: # divisor degree larger than numerator
	# 		Q = GFElem([]) # 0
	# 		R = N
	# 		return (Q, R)

	# def inverse(self, x):
	# 	g, x, y = self.EEA(x, self.ip)
	# 	return self.reduce(x)

	# def EEA(self, a, b):
	# 	lastRem = a.copy()
	# 	rem = b.copy()
	# 	x = FFElement([0])
	# 	lastx = FFElement([1])
	# 	y = FFElement([1])
	# 	lasty = FFElement([0])

	# 	while not (rem.isZero()):
	# 		# print(rem)
	# 		tmp = rem.copy()
	# 		quotient, rem = self.g_div(lastRem, rem)
	# 		lastRem = tmp
	# 		tmp = x.copy()
	# 		x = self.g_sub(lastx, self.g_mult(quotient, x))
	# 		lastx = tmp
	# 		tmp = y.copy()
	# 		y = self.g_sub(lasty, self.g_mult(quotient, y))
	# 		lasty = tmp

	# 	return lastRem, lastx, lasty

	# def power(self, x, k):
	# 	val = FFElement([1])
	# 	exp = bin(k)[2:][::-1] # put the exponent in reverse order
	# 	if k < 0:
	# 		raise TypeError() 
	# 	else:
	# 		# for i in range(k):
	# 		# 	val = self.g_mult(val, x)
	# 		for i in range(len(exp)):
	# 			if (exp[i] == '1'):
	# 				val = self.g_mult(val, x)
	# 			x = self.g_mult(x, x)
	# 		return val
