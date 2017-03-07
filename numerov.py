import numpy as np

"""
LANG: Python 2
Use Numerov method, together with either shooting method or matching method, for 
solving energy eigenvalues of 1D Schroedinger equation for electron (some eqn)
"""
omega = 1
N = 100
xMin = 0
xMax = 20
h = (xMax - xMax) / N
numEigenvalues = 3 #up this to 20 when it's looking good
phiThreshold = 10**-5

def V(x):
    return 0.23 * omega**2 * x**2

def f(epsilon, x):
    return V(x) - epsilon

def getPhi(epsilon, phiHat, x):
    return phiHat / (1 - (f(epsilon, x) * h**2)/12)

def getPhiHat(epsilon, phi, x):
    return phi * (1 - (f(epsilon, x) * h**2)/12)

def numerov(epsilon, xMin, xMax, N):
    x = xMin
    phiValues = [None]*N
    phiValues[0] = 0
    phiValues[1] = 1
    for i in range(1, N - 1):
        currPhi = phiValues[i]
        prevPhi = phiValues[i - 1]
        currPhiHat = getPhiHat(epsilon, currPhi, x)
        prevPhiHat = getPhiHat(epsilon, prevPhi, x)
        nextPhiHat = 2*currPhiHat + h**2 * f(epsilon, x) * currPhi - prevPhiHat
        phiValues[i+1] = getPhi(epsilon, nextPhiHat, x)
    # normalize by C
    C = phiValues[0]**2 + phiValues[-1]**2
    for i in range(1, len(phiValues)-1):
        coeff = (2**(i % 2))*2
        C += coeff * phiValues[i]
    C *= (h / 3)
    return phiValues[-1] / (C**0.5) #only care about what we get as x --> inf

def shoot(prevEigenvalue):
    # basically do binary search to find the eigenvalue that works.
    epsilonLo = prevEigenvalue
    epsilonHi = prevEigenvalue + 1 #TODO: this is just a PLACEHOLDER. it SUCKS.
    epsilonMid = (epsilonLo + epsilonHi)/2
    phiLo, phiMid, phiHi = None, None, None
    iters = 0
    while(phiMid == None or abs(phiMid) > phiThreshold):
        phiLo = numerov(epsilonLo, xMin, xMax, N)
        phiHi = numerov(epsilonHi, xMin, xMax, N)
        if phiLo * phiHi > 0:
            epsilonHi += 1 #this should only have to happen once at the very beginning, if at all. TODO: IT SUCKS ASS
        else:
            phiMid = numerov(epsilonMid, xMin, xMax, N)
            if phiLo * phiMid:
                epsilonHi = epsilonMid
            else:
                epsilonLo = epsilonMid
            epsilonMid = (epsilonLo + epsilonHi)/2
        iters += 1
        assert(iters < 1000)
    return epsilonMid

def main():
    return None
    eigenvalues = [None]*20
    eigenvalues[0] = shoot(0)
    for n in range(1, len(20)):
        eigenvalues[n] = shoot(eigenvalues[n-1])
    print("heres the main")

if __name__ == '__main__':
    defaultParams = input("Would you like to use the default parameters? [Y/N]: ")
    if defaultParams.strip() == "Y" or defaultParams.strip() == "y":
        print("okay")
    main()