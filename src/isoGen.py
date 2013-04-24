# File: isoGen.py
# Author: Christopher Wood

import sys
import pickle
from numpy import matrix # for matrix manipulation
sys.path.append("./cfa/")
from galois import GFElem, GFExtensionElem, GF, GFExtension

# A package containing the isomorphic mapping data so it can be serialized to a file
class IsoPackage:
    def __init__(self, neField, neIp, eField, eIp, alpha, beta, order, isoMap, T, Ti):
        self.neField = neField
        self.neIp = neIp
        self.eField = eField
        self.eIp = eIp
        self.alpha = alpha
        self.beta = beta
        self.order = order
        self.isoMap = isoMap
        self.T = T
        self.Ti = Ti

    def cost(self):
        gates = 0
        for i in range(self.neField.getExp()):
            for j in range(self.neField.getExp()):
                if (int(self.T[(i,j)]) == 1):
                    gates = gates + 1
        for i in range(self.neField.getExp()):
            for j in range(self.neField.getExp()):
                if (int(self.Ti[(i,j)]) == 1):
                    gates = gates + 1
        return gates

    def getBundle(self):
        return {"neField" : self.neField, "neIp" : self.neIp, "eField" : self.eField, "eIp" : self.eIp, "alpha" : self.alpha, "beta" : self.beta, "order" : self.order, "isoMap" : self.isoMap, "T" : self.T, "Ti" : self.Ti}
 
# Wrapper to hold an isomorphic map and easily display its contents
class IsoMap:
    def __init__(self, isoMap):
        self.isoMap = isoMap

    def __str__(self):
        result = ""
        for k in self.isoMap.keys():
            result = result + str(k) + " -> " + str(self.isoMap[k]) + "\n"
        return result

# Exhaustive search for isomorphic function being generated
def genIsomorphicFunction(neField, neIp, eField, eIp, alpha, beta, order):
    # alpha = GFElem
    # beta = GFExtensionElem
    isoMap = {}
    alphaMap = {}
    betaMap = {}
    for p in range(order):
        c1 = GFElem([]) # 0x^0
        c2 = GFExtensionElem([GFElem([])]) # 0y^0
        if (p != 0):
            c1 = alpha.copy()
            c2 = beta.copy()
            for i in range(1, p):
                c1 = neField.g_mult(c1, alpha) # multiply
                c2 = eField.g_mult(c2, beta) # multiply
        # c1 = c1 ^ i
        # c2 = c2 ^ i
        isoMap[str(c1)] = c2 # append the mapping (c1 will be unique if alpha is a generator!)
        alphaMap[str(c1)] = p # don't want to solve discrete log to find p :-)
        betaMap[str(c2)] = p # ditto

    # alpha^i -> beta^i preserves multiplicative isomorphism
    # now check to see if additive homomorphism is maintained
    firstOne = GFElem([1])
    secondOne = GFExtensionElem([GFElem([1])])
    for p in range(order):
        c1 = GFElem([]) # 0x^0
        c2 = GFExtensionElem([GFElem([])]) # 0y^0
        if (p != 0):
            c1 = alpha.copy()
            c2 = beta.copy()
            for i in range(1, p):
                c1 = neField.g_mult(c1, alpha) # multiply
                c2 = eField.g_mult(c2, beta) # multiply
        # c1 = c1 ^ i
        # c2 = c2 ^ i
        # now add 1 to both of these and then check the mapping...
        c1p1 = neField.g_add(c1, firstOne) # alpha^i + 1
        c2p1 = eField.g_add(c2, secondOne) # beta^i + 1

        # Pull r, the exponent of \alpha^r = \alpha^i + 1
        r1 = alphaMap[str(c1p1)]
        r2 = betaMap[str(c2p1)]
        if (r1 != r2): # if (alpha^r = alpha^i + 1)
            return None

    # Return the valid mapping
    return isoMap

def findMappingsFormal(neField, neIp, eField, eIp, neGens, eGens, outFile):
    order = neField.getOrder()
    mapIndex = 0
    maps = []
    for g1 in neGens:
        for g2 in eGens: 
            print >> sys.stderr, 'Checking mapping from ' + str(g1) + ' to ' + str(g2)
            isoMap = genIsomorphicFunction(neField, neIp, eField, eIp, g1, g2, order)
            if (isoMap != None):
                # Generate the transformation matrices using the two basis elements
                aPowers = []
                alpha = GFExtensionElem([GFElem([1])])
                for i in range(neField.getExp()):
                    aPowers.insert(0, (alpha, alpha.toBinary(eField.getBaseField().getExp(), eField.getExtension())))
                    alpha = eField.g_mult(alpha, g2)

                T = []
                for i in range(neField.getExp()):
                    T.append([])
                    for j in range(neField.getExp()):
                        T[i].append(aPowers[j].toBinary(eField.getN(), eField.getM())[i])
                Tm = matrix(T)
                Tim = Tm.I
                print >> sys.stderr, str(Tm)
                print >> sys.stderr, str(Tim)
                maps.append[(Tm, Tim)] # app

                # # Now build, store, and output the matrices.
                # T = {}
                # Ti = {} #Ti is just the INVERSE of T!
                # for i in range(neField.getExp()):
                #     for j in range(neField.getExp()):
                #         T[(i,j)] = aPowers[j][1][i]
                #         Ti[(i,j)] = bPowers[j][1][i]
                # print >> sys.stderr, "----"
                # for i in range(neField.getExp()):
                #     print >> sys.stderr, "[",
                #     for j in range(neField.getExp()):
                #         print >> sys.stderr, str(T[(i,j)]),
                #     print >> sys.stderr, "]"
                # print >> sys.stderr, "----"
                # for i in range(neField.getExp()):
                #     print >> sys.stderr, "[",
                #     for j in range(neField.getExp()):
                #         print >> sys.stderr, str(Ti[(i,j)]),
                #     print >> sys.stderr, "]"
                # print >> sys.stderr, "----"

                # Dump the package
                package = IsoPackage(neField, neIp, eField, eIp, g1, g2, order, isoMap, Tm, Tim)
                output = open(outFile + "_" + str(mapIndex) + '.pkl', 'wb') 
                mapIndex = mapIndex + 1
                pickle.dump(package, output)
                output.close()

                maps.append((mapWrap, T))
    return maps

def findMappingsExhaustive(neField, neIp, eField, eIp, neGens, eGens, outFile):
    order = neField.getOrder()
    mapIndex = 0
    maps = []
    for g1 in neGens:
        for g2 in eGens: 
            print >> sys.stderr, 'Checking mapping from ' + str(g1) + ' to ' + str(g2)
            isoMap = genIsomorphicFunction(neField, neIp, eField, eIp, g1, g2, order)
            if (isoMap != None):
                # Give some feedback...
                print >> sys.stderr, 'Found a map from ' + str(g1) + ' to ' + str(g2)

                # Display the map
                mapWrap = IsoMap(isoMap)
                print >> sys.stderr, str(mapWrap)

                basisElems = []
                for i in range(neField.getExp()):
                    elem = GFElem([], id = (1 << i), base = neField.getBase(), exp = neField.getExp())
                    basisElems.insert(0, isoMap[str(elem)])
                T = []
                for i in range(neField.getExp()):
                    T.append([])
                    for j in range(neField.getExp()):
                        T[i].append(basisElems[j].toBinary(eField.getN(), eField.getM())[i])
                Tm = matrix(T)
                Tim = Tm.I

                # Convert everything to positive integers (numpy throws in negatives)
                for i in range(neField.getExp()):
                    for j in range(neField.getExp()):
                        Tim[i,j] = int(abs(Tim[i,j]))

                # Dump the package
                package = IsoPackage(neField, neIp, eField, eIp, g1, g2, order, isoMap, Tm, Tim)
                output = open(outFile + "_" + str(mapIndex) + '.pkl', 'wb') 
                mapIndex = mapIndex + 1
                pickle.dump(package, output)
                output.close()

                print >> sys.stderr, str(Tm)
                print >> sys.stderr, str(Tim)
                maps.append((Tm, Tim, package.cost())) # app
    return maps

# Usage: python isoGen.py file
# Example file: 10011 2 4 111 2 2 010110 2

# Main entry point for test cases/isomorphic function generation
def main():
    # Open the file...
    f = open(sys.argv[1], "r")
    outFile = sys.argv[2]
    index = 0
    line = f.readline()
    while len(line) > 0:
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

        # Read in the generators for the base and extension field...
        neGensElems = f.readline()
        neGens = []
        data = neGensElems.split(" ")
        for g in data:
            g = g.strip()
            if (len(g) > 0):
                elem = GFElem([], binString = g)
                print(elem)
                neGens.append(elem)
        eGensElems = f.readline()
        eGens = []
        data = eGensElems.split(" ")
        for g in data:
            g = g.strip()
            if (len(g) > 0):
                elem = GFExtensionElem([], binString = g, smallSize = eExp)
                print(elem)
                eGens.append(elem)
        
        # neGens1 = neField.findGenerators()
        # eGens1 = eField.findGenerators()
        maps = findMappingsExhaustive(neField, neIp, eField, eIp, neGens, eGens, outFile + "_" + str(index))

        # Carry over to the next line in the file...
        line = f.readline()
        index = index + 1
	
if __name__ == '__main__':
	main()