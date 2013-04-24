# File: inverseComplexity.py
# Author: Christopher Wood

from galois import *

def complexity(A, B):
	mults = 0
	weight = 0
	if not A.isUnit():
		mults = mults + 2
		weight = weight + A.wt()
	if not B.isUnit():
		mults = mults + 1
		weight = weight + B.wt()
	return mults, weight

def main():
	print(complexity(GFElem([0,0,0,1]), GFElem([0,0,1,0])))

if __name__ == "__main__":
	main()