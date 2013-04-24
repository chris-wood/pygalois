# File: cfa.py
# Author: Christopher Wood, caw4567@rit.edu

import sys
import pickle
from galois import GFElem, GFExtensionElem, GF, GFExtension

# A package containing the isomorphic mapping data so it can be serialized to a file
class IsoPackage:
    def __init__(self, neField, neIp, eField, eIp, alpha, beta, order, isoMap):
        self.neField = neField
        self.neIp = neIp
        self.eField = eField
        self.eIp = eIp
        self.alpha = alpha
        self.beta = beta
        self.order = order
        self.isoMap = isoMap

    def getBundle(self):
        return {"neField" : self.neField, "neIp" : self.neIp, "eField" : self.eField, "eIp" : self.eIp, "alpha" : self.alpha, "beta" : self.beta, "order" : self.order, "isoMap" : self.isoMap}
 
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
            c1 = alpha
            c2 = beta
            for i in range(1, p):
                c1 = neField.g_mult(c1, alpha) # multiply
                c2 = eField.g_mult(c2, beta) # multiply
        # c1 = c1 ^ i
        # c2 = c2 ^ i
        #isoMap.append((c1, c2)) # append the mapping
        isoMap[str(c1)] = c2 # append the mapping (c1 will be unique if alpha is a generator!)
        alphaMap[str(c1)] = p # don't want to solve discrete log to find p :-)
        betaMap[str(c2)] = p # don't want to solve discrete log to find p :-)

    # debug!
    #print(isoMap)

    # alpha^i -> beta^i preserves multiplicative isomorphism
    # now check to see if additive homomorphism is maintained
    firstOne = GFElem([1])
    secondOne = GFExtensionElem([GFElem([1])])
    for p in range(order):
        c1 = GFElem([]) # 0x^0
        c2 = GFExtensionElem([GFElem([])]) # 0y^0
        if (p != 0):
            c1 = alpha
            c2 = beta
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

        # if not (isoMap[str(c1p1)] == c2p1): # if the additive group homomorphism doesn't stick, return an invalid mapping
        #     return None # return None as invalid mapping

    # Return the valid mapping
    return isoMap

def genSmallCase():
    neIp = GFElem([1,0,0,1,1]) # x^4 + x + 1
    neField = GF(2, 4, neIp)
    print("BIG FIELD")
    print(neField)

    print("LISTING ALL OF THE BIG FIELD ORDERS")
    gens1 = []
    order = (2**4) - 1
    # print(field.findOrder(0x03)[1])
    for i in range(order + 1):
        result = neField.findOrder(i)
        print(i, result[1])
        if (result[1] == order):
            gens1.append(result[0])
    # print("LISTING ALL OF THE GENERATORS (as polynomials)")
    print("(" + str(len(gens1)) + " generators total)")
    for g in gens1:
        print(g)

    ip = GFElem([1,1,1]) # x^2 + x + 1
    smallField = GF(2, 2, ip)
    eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,0])]) # x^2 + x + 10
    eField = GFExtension(smallField, 2, eIp)
    print("SMALL FIELD")
    print(str(eField))

    print("LISTING ALL OF THE SMALL FIELD ORDERS")
    gensExt = []
    for i in range(order + 1):
        result = eField.findOrder(i)
        if (result[1] == (order)): # it's a generator...
            gensExt.append(result[0]) 
        print(str(result[0]) + ": order = " + str(result[1]))
    print("Total number of generators: " + str(len(gensExt)))
    for alpha in gensExt:
        print(str(alpha))

    print("FINDING ALL POSSIBLE ISOMORPHIC MAPPINGS")
    order = 2 ** 4
    characteristic = 2
    exponent = 4
    subExponent = 2
    mapIndex = 0
    for g1 in gens1:
        for g2 in gensExt: # TODO: BRING BACK WHEN TECHNIQUE IS DEEMED CORRECT
            print("Checking mapping from " + str(g1) + " to " + str(g2))
            isoMap = genIsomorphicFunction(neField, neIp, eField, eIp, g1, g2, order)
            if (isoMap != None):
                print("FOUND IT")
                # output = open('smallAlpha.pkl', 'wb') #alpha.pkl
                # pickle.dump(g1, output)
                # output.close()
                # output = open('smallBeta.pkl', 'wb') #beta.pkl
                # pickle.dump(g2, output)
                # output.close()
                # output = open('smallIsomap.pkl', 'wb')
                # pickle.dump(isoMap, output)
                # output.close()

                # Dump the package
                package = IsoPackage(neField, neIp, eField, eIp, g1, g2, order, isoMap)
                output = open('smallBundle_isoMap_' + str(mapIndex) + '.pkl', 'wb') 
                mapIndex = mapIndex + 1
                pickle.dump(package, output)
                output.close()

                # Display the map
                mapWrap = IsoMap(isoMap)
                print(mapWrap)

                # Generate the transformation matrix
                powers = []
                alpha = GFExtensionElem([GFElem([1])])
                for i in range(neField.getExp()):
                    powers.insert(0, (alpha, alpha.toBinary(eField.getBaseField().getExp(), eField.getExtension())))
                    alpha = eField.g_mult(alpha, g2)
                T = {}
                for i in range(neField.getExp()):
                    for j in range(neField.getExp()):
                        T[(i,j)] = powers[j][1][i]
                for i in range(neField.getExp()):
                    print("["),
                    for j in range(neField.getExp()):
                        print(str(T[(i,j)])),
                    print("]\n")
                return T

def findMappings(neField, neIp, eField, eIp, neGens, eGens):
    order = neField.getOrder()
    mapIndex = 0
    maps = []
    for g1 in neGens:
        # n = eField.getBaseField().getExp()
        # m = eField.getExtension()
        # r = ((2 ** (n * m)) - 1) / ((2 ** n) - 1)
        # g
        # for r = 0
        for g2 in eGens: 
            #g2 = GFExtensionElem([GFElem([1,0]), GFElem([1,1,1,0])])
            print >> sys.stderr, 'Checking mapping from ' + str(g1) + ' to ' + str(g2)
            isoMap = genIsomorphicFunction(neField, neIp, eField, eIp, g1, g2, order)
            if (isoMap != None):
                # Dump the package
                package = IsoPackage(neField, neIp, eField, eIp, g1, g2, order, isoMap)
                output = open('isoMap_' + str(mapIndex) + '.pkl', 'wb') 
                mapIndex = mapIndex + 1
                pickle.dump(package, output)
                output.close()

                # Display the map
                mapWrap = IsoMap(isoMap)
                print("****FOUND ONE MAP****")
                print(mapWrap)

                # Generate the transformation matrix
                powers = []
                alpha = GFExtensionElem([GFElem([1])])
                for i in range(neField.getExp()):
                    powers.insert(0, (alpha, alpha.toBinary(eField.getBaseField().getExp(), eField.getExtension())))
                    alpha = eField.g_mult(alpha, g2)
                T = {}
                for i in range(neField.getExp()):
                    for j in range(neField.getExp()):
                        T[(i,j)] = powers[j][1][i]
                for i in range(neField.getExp()):
                    print("["),
                    for j in range(neField.getExp()):
                        print(str(T[(i,j)])),
                    print("]\n")
                print("********")
                maps.append((mapWrap, T))
                return maps
    return maps

def findGenerators(field):
    fieldOrder = field.getOrder()
    gens = []
    for i in range(fieldOrder):
        elemOrder = field.findOrder(i)
        if (elemOrder[1] == (fieldOrder - 1)): # it's a generator...
            gens.append(elemOrder[0]) 
    return gens

# Main entry point for test cases/isomorphic function generation
def main():
    genSmallCase()

    return # quick test above

    # COMMENTED OUT FOR DEMO
    # Rijndael GF - GF(2^8) defined by the IP below
    # neIp = GFElem([1,0,0,0,1,1,1,0,1]) #x^8 + x^4 + x^3 + x^2 + 1 # CHANGED THIS TO PRIMITIVE POLYNOMIAL
    # neField = GF(2, 8, neIp) # GF(2^8)
    # ip = GFElem([1,0,0,1,1]) # x^4 + x + 1
    # smallField = GF(2, 4, ip)
    # eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,1,0,0])]) # x^2 + x + 1100
    # eField = GFExtension(smallField, 2, eIp) # GF((2^4)^2)

    # Small field test case
    # neIp = GFElem([1,0,0,1,1]) # x^4 + x + 1
    # neField = GF(2, 4, neIp)
    # ip = GFElem([1,1,1]) # x^2 + x + 1
    # smallField = GF(2, 2, ip)
    # eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,0])]) # x^2 + x + 10
    # eField = GFExtension(smallField, 2, eIp)

    # Find all of the generators for the fields
    print("Finding generators")
    neGens = neField.findGenerators()
    eGens = eField.findGenerators()

    # Search for all isomorphic mappings
    print(neGens)
    print(eGens)
    print("Finding the isomorphic mappings")
    maps = findMappings(neField, neIp, eField, eIp, neGens, eGens)

# Some test cases for the programs
def test():
    ip = GFElem([1,0,0,0,1,1,0,1,1]) #x^8 + x^4 + x^3 + x + 1
    field = GF(2, 8, ip)

    # Sample polynomials
    x = GFElem([1,1])
    y1 = GFElem([1,1])
    y2 = GFElem([1,0])
    y3 = GFElem([])
    y4 = GFElem([1,0,0,0,0,0,0,0])
    y5 = GFElem([1,0,0,0,0,0,0])
    v3 = GFElem([1,1])
    v4 = GFElem([1,0])
    v7 = GFElem([1,1,1])
    x10 = GFElem([1,1,1,1,1,1,1,1])
    y10 = GFElem([0,0,0,0,0,0,1,1])
    x11 = GFElem([0,1,0,1,0,0,1,1])
    y11 = GFElem([1,1,0,0,1,0,1,0])
    
    print("STARTING TEST CASES FOR GF(2^8)")

    # Sum test cases
    sum = field.g_add(x, y1)
    print(str(x) + " + " + str(y1) + " = " + str(sum))
    sum = field.g_add(x, y2)
    print(str(x) + " + " + str(y2) + " = " + str(sum))

    # Diff test cases
    diff = field.g_sub(x, y1)
    print(str(x) + " - " + str(y1) + " = " + str(diff))
    diff = field.g_sub(x, y2)
    print(str(x) + " - " + str(y2) + " = " + str(diff))

    # Mult test cases
    mult = field.g_mult(x, y1)
    print(str(x) + " * " + str(y1) + " = " + str(mult))
    mult = field.g_mult(v3, v7)
    print(str(v3) + " * " + str(v7) + " = " + str(mult))
    mult = field.g_mult(y5, y5)
    print(str(y5) + " * " + str(y5) + " = " + str(mult))
    mult = field.g_mult(y4, y4)
    print(str(y4) + " * " + str(y4) + " = " + str(mult))
    mult = field.g_mult(x10, y10)
    print(str(x10) + " * " + str(y10) + " = " + str(mult))
    mult = field.g_mult(x11, y11)
    print(str(x11) + " * " + str(y11) + " = " + str(mult))

    # Division cases
    div = field.g_div(x, ip)
    print(str(x) + " / " + str(ip) + " = " + str(div[0]) + "," + str(div[1]))

    print("STARTING TEST CASES FOR GF(2^2)")

    ip = GFElem([1, 1, 1]) #x^2 + x + 1
    field = GF(2, 2, ip)
    print(field)
    x0 = GFElem([0, 0])
    x1 = GFElem([0, 1])
    x2 = GFElem([1, 0])
    x3 = GFElem([1, 1])
    mult = field.g_mult(x0, x1)
    print(str(x0) + " * " + str(x1) + " = " + str(mult))
    mult = field.g_mult(x0, x2)
    print(str(x0) + " * " + str(x2) + " = " + str(mult))
    mult = field.g_mult(x0, x3)
    print(str(x0) + " * " + str(x3) + " = " + str(mult))
    mult = field.g_mult(x1, x2)
    print(str(x1) + " * " + str(x2) + " = " + str(mult))
    mult = field.g_mult(x1, x3)
    print(str(x1) + " * " + str(x3) + " = " + str(mult))
    mult = field.g_mult(x2, x3)
    print(str(x2) + " * " + str(x3) + " = " + str(mult))

    print("STARTING TEST CASES FOR GF(4^2)")

    ip = GFElem([1, 1, 1]) #x^2 + x + 1
    field = GF(4, 2, ip)
    print(field)
    x0 = GFElem([0, 0])
    x1 = GFElem([0, 1])
    x2 = GFElem([1, 0])
    x3 = GFElem([1, 1])
    mult = field.g_mult(x0, x1)
    print(str(x0) + " * " + str(x1) + " = " + str(mult))
    mult = field.g_mult(x0, x2)
    print(str(x0) + " * " + str(x2) + " = " + str(mult))
    mult = field.g_mult(x0, x3)
    print(str(x0) + " * " + str(x3) + " = " + str(mult))
    mult = field.g_mult(x1, x2)
    print(str(x1) + " * " + str(x2) + " = " + str(mult))
    mult = field.g_mult(x1, x3)
    print(str(x1) + " * " + str(x3) + " = " + str(mult))
    mult = field.g_mult(x2, x3)
    print(str(x2) + " * " + str(x3) + " = " + str(mult))
    
if (__name__ == "__main__"):
    main()
