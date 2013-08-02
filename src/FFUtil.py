# TODO: implement Miller-Rabin and AKS primality tests

def slowFactor(n):
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
			return d 
	return n