# File: polyGen.py
# Author: Christopher Wood
# Description: Generates irreducible and primitive polynomials

import sys
from galois import *
from galois_util import *

def findIrreduciblePolynomialsBase(base, n, lb, ub):
	order = base ** n
	coeff = []
	irrPolys = []
	primPolys = []
	for i in range(n + 1):
		coeff.append(0)
	# fieldOrder = base ** (n+1)
	for i in range(lb, ub + 1):
		coeff = createBaseElem(base, i, n + 1)
		elem = GFElem(coeff)
		f = GF(base, n, elem)
		if (elem.degree() == n and elem[elem.degree()] == 1 and elem[0] == 1): # monic polynomials only
			print >> sys.stderr, "Index: " + str(i)
			print >> sys.stderr, "Trying: " + str(elem)
			try:
				gens = f.findGenerators(earlyTerm = True)
				if (len(gens) > 0):
					print >> sys.stderr, "Irreducible"
					irrPolys.append(elem)
					x = GFElem([1,0])
					if (x in gens):
						print >> sys.stderr, "Primitive"
						primPolys.append(elem)
			except Exception as e:
				print >> sys.stderr, str(e)
				pass
	return irrPolys, primPolys

def findIrreduciblePolynomialsExtension(base, n, m, lb, ub, simple = False):
	order = base ** (n + m)
	coeff = []
	bps = []
	basePolys, throwaway = findIrreduciblePolynomialsBase(base, n, 0, (base ** (n + 1)) - 1)
	print >> sys.stderr, "Trying all combinations..."
	index = 0
	result = []
	for bp in basePolys:
		print >> sys.stderr, "Base polynomial: " + str(bp.getCoeff())
		
		# Lists to store each set of polynomials here...
		irrPolys = []
		primPolys = []
		
		# Create the base field
		ip = GFElem(bp.getCoeff())
		baseField = GF(2, n, ip)
		coeff = []
		for i in range(m  + 1):
			coeff.append(GFElem([0]))
		fieldOrder = (base ** n) ** m # m is the extension
		elems = []
		
		# Now try every polynomial in the larger field
		for i in range(lb, ub + 1): # caw: had this as -1 before
			p = []
			subOrder = baseField.getOrder()
			for j in range(m):
				gg = (i >> (j * n))
				bits = []
				for b in range(n):
					bits.append((gg & (1 << (n - b - 1))) >> (n - b - 1))
				coeff[j+1] = GFElem(bits)
			coeff[0] = GFElem([1]) # only check minimal polynomials

			# Make sure all terms except the last are set to one (special case of polynomial generation)
			if (simple):
				for i in range(1, m):
					coeff[i] = GFElem([1]) 

			# Create the candidate polynomial and then check its properties before proceeding
			elem = GFExtensionElem(coeff)
			if not (elem.bin(n) in elems): # don't visit a polynomial we've already seen...
				if (not elem[0].isZero() and not elem[len(elem.getCoeff()) - 1].isZero()):
					elems.append(elem.bin(n))
					f = GFExtension(baseField, m, elem)
					print >> sys.stderr, "Index: " + str(i)
					print >> sys.stderr, "Trying: " + str(elem)
					try:
						gens = f.findGenerators(earlyTerm = True)
						if (len(gens) > 0):
							c = coeff[:]
							if (c not in irrPolys):
								irrPolys.append(c)
							print >> sys.stderr, "Irreducible."
							x = GFExtensionElem([GFElem([1]), GFElem([0])])
							if (x in gens):
								if (c not in primPolys):
									print >> sys.stderr, "Primitive."
									primPolys.append(c)
					except Exception as e:
						print >> sys.stderr, "Error in findIrreduciblePolynomialsExtension: " + str(e)
						pass
				else:
					# print >> sys.stderr, "Polynomial didn't pass the checkpoint"
					pass
			else:
				# print >> sys.stderr, "Duplicate detected."
				pass
		# Make copies of irrPolys and primPolys and then append to the output list
		result.append((bp, irrPolys[:], primPolys[:]))
	return result

def smallTest():
	print("-----------------------------------------------")
	print("For field GF(2^4), irreducible polynomials are:")
	irrPolys, primPolys = findIrreduciblePolynomialsBase(2, 4, 0, 31)
	for p in irrPolys:
		elem = GFElem(p)
		print(str(elem.bin()))
	print("Primitive polynomials are:")
	for p in primPolys:
		elem = GFElem(p)
		print(str(elem.bin()))
	print("-----------------------------------------------")
	print("For field GF((2^2)^2), irreducible polynomials are:")
	basePolys, irrPolys, primPolys = findIrreduciblePolynomialsExtension(2, 2, 2, 0, 31)
	index = 0
	for bp in basePolys:
		print("Base polynomial:")
		baseElem = GFElem(bp)
		print(baseElem)
		print("Irreducibles:")
		for p in irrPolys[index]:
			elem = GFExtensionElem(p)
			print(elem.bin(2))
		index = index + 1
	print("Primitive polynomials are:")
	index = 0
	for bp in basePolys:
		print("Base polynomial:")
		baseElem = GFElem(bp)
		print(baseElem)
		print("Primitives:")
		for p in primPolys[index]:
			elem = GFExtensionElem(p)
			print(str(elem.bin(2)))
		index = index + 1

def aesTest():
	print("-----------------------------------------------")
	print("For field GF((2^4)^2), irreducible polynomials are:")
	result = findIrreduciblePolynomialsExtension(2, 4, 2, 0, 511, simple = True)
	index = 0
	n = 4
	for r in result:
		bp = r[0]
		ips = r[1]
		ps = r[2]
		print("ip : " + str(bp.bin()) + " : "),
		for p in ips:
			elem = GFExtensionElem(p)
			print(elem.bin(n) + ","),
		print("p : " + str(bp.bin()) + " : "),
		for p in ps:
			elem = GFExtensionElem(p)
			print(elem.bin(n) + ","),
		print("")

def specificTest():
	ip = GFElem([1,0,0,1,1]) # x^4 + x + 1
	smallField = GF(2, 4, ip)
	eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,0,0,1])]) # x^2 + x + 1100
	eField = GFExtension(smallField, 2, eIp) # GF((2^4)^2)
	gens = eField.findGenerators()
	print(gens)

def test():
	aesTest()

# Usage: python polyGen.py MODE (options)
# MODE = 0, base n lb ub
# MODE = 1, base n m lb ub
# MODE = 2, run test method

def main():
	ext = int(sys.argv[1])
	if (ext == 0):
		base = int(sys.argv[2])
		n = int(sys.argv[3])
		lb = int(sys.argv[4])
		ub = int(sys.argv[5])
		print >> sys.stderr, "GF(" + str(base) + "^" + str(n) + ")"
		irrPolys, primPolys = findIrreduciblePolynomialsBase(base, n, lb, ub)
		print(len(irrPolys))
		for p in irrPolys:
			elem = GFElem(p)
			print(elem.bin() + ","),
		print("")
		print(len(primPolys))
		for p in primPolys:
			elem = GFElem(p)
			print(elem.bin() + ","),
		print("")
	elif (ext == 1):
		base = int(sys.argv[2])
		n = int(sys.argv[3])
		m = int(sys.argv[4])
		lb = int(sys.argv[5])
		ub = int(sys.argv[6])
		print >> sys.stderr, "GF((" + str(base) + "^" + str(n) + ")^" + str(m) + ")"
		result = findIrreduciblePolynomialsExtension(base, n, m, lb, ub)
		for r in result:
			bp = r[0]
			ips = r[1]
			ps = r[2]
			print("ip : " + str(bp.bin()) + " : "),
			for p in ips:
				elem = GFExtensionElem(p)
				print(elem.bin(n) + ","),
			print("p : " + str(bp.bin()) + " : "),
			for p in ps:
				elem = GFExtensionElem(p)
				print(elem.bin(n) + ","),
			print("")
		# index = 0
		# for bp in basePolys:
			# print(bp.bin() + ","),
			#print("Irreducibles:")
			# for p in irrPolys[index]:
				# elem = GFExtensionElem(p)
				# print(elem.bin(n) + ","),
			# print("")
			# index = index + 1
		# index = 0
		# for bp in basePolys:
			# print(bp.bin() + ","),
			#print("Irreducibles:")
			# for p in primPolys[index]:
				# elem = GFExtensionElem(p)
				# print(elem.bin(n) + ","),
			# print("")
			# index = index + 1
	elif (ext == 2):
		base = int(sys.argv[2])
		k = int(sys.argv[3])
		n = int(sys.argv[4])
		m = int(sys.argv[5])
		lb = int(sys.argv[6])
		ub = int(sys.argv[7])
		print >> sys.stderr, "GF(" + str(base) + "^" + str(k) + ")"
		print >> sys.stderr, "GF((" + str(base) + "^" + str(n) + ")^" + str(m) + ")"
		irrPolys, primPolys = findIrreduciblePolynomialsBase(base, k, lb, ub)
		result = findIrreduciblePolynomialsExtension(base, n, m, lb, ub)
		for ip in irrPolys:
			for r in result:
				bp = r[0]
				ips = r[1]
				ps = r[2]
				for p in ips:
					elem = GFExtensionElem(p)
					print(ip.bin() + " " + str(base) + " " + str(k) + " " + bp.bin() + " " + str(n) + " " + elem.bin(n) + " " + str(m))
	elif (ext == 3):
		test()

if __name__ == "__main__":
	main()