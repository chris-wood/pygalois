# File: galois_util.py
# Author: Christopher Wood

from galois import *

def createBaseElem(coeff, base, index):
	ii = index
	for j in range(len(coeff)):
		coeff[j] = ii % base
		ii = ii / base
	elem = GFElem(coeff)
	return elem