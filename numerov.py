import numpy as np


# LANG: Python 3
# Use Numerov method, together with either shooting method or matching method, for 
# solving energy eigenvalues of 1D Schroedinger equation for electron (some eqn)

def externalPotential(x, omega=1):
    return 0.25 * (omega**2) * (x**2)

def f(epsilon, x, omega=1):
    return externalPotential(x, omega=omega) - epsilon

def numerov(epsilon, omega=1, xMax = 1, xMin=0, N=1000):
    h = abs(xMax - xMin)/N
    xVals = [xMin + i*h for i in range(N+1)]
    fVals = [f(epsilon, xVals[i], omega) for i in range(N+1)]
    phiVals = [None]*(N+1)
    phiVals[0] = 0
    phiVals[1] = 1
    for i in range(1, N):
        prevPhiHat = (1 - (1/12)*(h**2)*fVals[i-1])*phiVals[i-1]
        nextPhiHat = (2 + (5/6)*(h**2)*fVals[i])*phiVals[i] - prevPhiHat
        # nextPhiHat = (2 - (5/6)*(h**2)*fVals[i])*phiVals[i] - prevPhiHat #idk, i'm just pulling this from some pdf
        phiVals[i+1] = nextPhiHat / (1 - (1/12)*(h**2)*fVals[i+1])
    return phiVals

def getNormalizedPhiVals(phiVals, xMax, xMin, N):
    h = abs(xMax - xMin)/N
    C = 0
    for i in range(1, N):
        C += (2*h/3)*(1 + (i%2)) * phiVals[i]**2
    C += (h/3)*phiVals[0] + (h/3)*phiVals[N]
    for i in range(len(phiVals)):
        phiVals[i] = phiVals[i]/np.sqrt(C)
    return phiVals

def getPhiVals(epsilon, omega=1, xMax = 1, xMin=0, N=1000):
    unnormalizedPhiVals = numerov(epsilon, omega=omega, xMax=xMax, xMin=xMin, N=N)
    return getNormalizedPhiVals(unnormalizedPhiVals, xMax=xMax, xMin=xMin, N=N)


def getWaveFunctionVals():
    numEpsilonVals = 25
    maxEpsilonVal = 5
    minEpsilonVal = 0
    epsilonVals = [(i+1)*(maxEpsilonVal - minEpsilonVal)/numEpsilonVals for i in range(numEpsilonVals)]
    
    xMin = -5
    xMax = 5
    N = 10000
    #defining the following 2 because useful for plotting
    h = abs(xMax - xMin)/N
    xVals = [xMin + i*h for i in range(N+1)]
    # 2D list. rows correspond to series from numerov recursion. columns 
    # correspond to column'th value in the recursion  over all epsilons
    phiVals = [None]*numEpsilonVals
    lastPhiVals = [None]*numEpsilonVals
    for i in range(len(epsilonVals)):
        phiVals[i] = getPhiVals(epsilonVals[i], xMin=xMin, xMax=xMax, N=N)
        lastPhiVals[i] = phiVals[i][-1]
    for i in range(len(lastPhiVals)):
        print(epsilonVals[i], lastPhiVals[i])
    plotWaveFunctionVals(xVals, epsilonVals, phiVals)

def plotWaveFunctionVals(xVals, epsilonVals, phiVals):
    try:
        import matplotlib.pyplot as plt 
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        X, Y = np.meshgrid(xVals, epsilonVals)
        Z = np.array(phiVals).reshape(X.shape)
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)
        ax.set_xlabel("x values")
        ax.set_ylabel("$\epsilon$ values")
        ax.set_zlabel("wave amplitudes")
        ax.set_zlim(0, 1)
        fig.colorbar(surf, shrink=0.5, aspect=10)
        plt.show()
    except:
        print("that's all folks")

def shoot(epsilonLo, epsilonHi, omega, xMin, xMax, N, shootingTolerance=10**-3):
    pass

def getEigenvaluesByShooting(omega=1, xMax = 1, xMin=0, N=1000):
    numEigenvalues = 2
    # Have guesses of eigenvalues - epsilon. 
    # Solve for values of epsilon such that phi(xMax; epsilon) close to 0.
    # Make no assumptions about eigenvalues except that they're all positive,
    # so incrementally increase epsilon until the sign of phi(xMax, epsilon)
    # changes. then use shoot with args epsilon and epsilon - epsilonStep to 
    # solve for eigenvalue. Then continue incrementing as if nothing happened.
    # Know eigenvalues will be positive.
    epsilon = 0 
    #sort of like the nyquist frequency
    epsilonStep = omega/10 
    eigenvalues = []
    prevPhiVal = numerov(epsilon, omega=omega, xMax = xMax, xMin=xMin, N=N)
    currPhiVal = None
    while(len(eigenvalues) < numEigenvalues):
        epsilon += epsilonStep
        currPhiVal = numerov(epsilon, omega=omega, xMax = xMax, xMin=xMin, N=N)
        # if we've encountered a sign change:
        # equivalent to currPhiVal * prevPhiVal < 0, but without overflow worry
        if ((currPhiVal > 0 and prevPhiHat < 0) or
                    (currPhiVal < 0 and prevPhiHat > 0)):
            shot = shoot(epsilonHi=epsilon, epsilonLo=(epsilon-epsilonStep), 
                omega=omega, xMax = xMax, xMin=xMin, N=N, 
                shootingTolerance=10**-1)
            eigenvalues.append(shot)


if __name__ == '__main__':
    getWaveFunctionVals()
    # pass