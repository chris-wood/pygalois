# File: genGen.py
# Author: Christopher Wood

import sys
import pickle
from numpy import matrix # for matrix manipulation
sys.path.append("./cfa/")
from galois import GFElem, GFExtensionElem, GF, GFExtension

# Main entry point for test cases/isomorphic function generation
def main():
	# Open the file...
	f = open(sys.argv[1], "r")
	for line in f:
		data = line.split(" ")

		# Read all of the field parameters...
		rz = str(data[0])
		neBase = int(data[1])
		neExp = int(data[2])
		px = str(data[3])
		eBase = int(data[4])
		eExp = int(data[5])
		qy = str(data[6])
		ext = int(data[7])

		neIp = GFElem([], binString = rz)
		neField = GF(neBase, neExp, neIp)
		smallIp = GFElem([], binString = px)
		smallField = GF(eBase, eExp, smallIp)
		eIp = GFExtensionElem([], binString = qy, smallSize = eExp)
		eField = GFExtension(smallField, ext, eIp)

		# Compute all of the generators...
		neGens = neField.findGenerators()
		eGens = eField.findGenerators()

		print(line),
		for g in neGens:
			print >> sys.stderr, str(g)
			print >> sys.stderr, str(g.bin())
			print(g.bin()),
		print("")
		for g in eGens:
			print >> sys.stderr, str(g)
			print >> sys.stderr, str(g.bin(eExp))
			print(g.bin(eExp)),
		print("")

# Entry point
if __name__ == '__main__':
	main()