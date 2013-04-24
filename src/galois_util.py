# File: galois_util.py
# Author: Christopher Wood

import math

def createBaseElem(base, index, dim):
	ii = index
	coeff = []
	for j in range(dim):
		coeff.append(ii % base)
		ii = ii / base
	coeff.reverse()
	return coeff

def primeFactors(n):
	factors = []
	while (n > 1):
		factor = getSinglePrimeFactor(n)
		if not (factor in factors): 
			factors.append(factor)
		n = n / factor
	return tuple(factors)

def getSinglePrimeFactor(n):
	if (n % 2 == 0):
		return 2
	for d in range(3, int(math.ceil(math.sqrt(n)) + 1), 2):
		if (n % d == 0):
			return d # divisor.
	return n

def main():
	print(primeFactors(255))

if __name__ == "__main__":
	main()