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
    phiVals[1] = h
    for i in range(1, N):
        prevPhiHat = (1 - (1/12)*(h**2)*fVals[i-1])*phiVals[i-1]
        nextPhiHat = (2 + (5/6)*(h**2)*fVals[i])*phiVals[i] - prevPhiHat
        # nextPhiHat = (2 - (5/6)*(h**2)*fVals[i])*phiVals[i] - prevPhiHat #idk, i'm just pulling this from some pdf
        phiVals[i+1] = nextPhiHat / (1 - (1/12)*(h**2)*fVals[i+1])
    return phiVals

def normalizedPhiVals(phiVals, xMax, xMin, N):
    h = abs(xMax - xMin)/N
    C = 0
    for i in range(1, N):
        C += (2*h/3)*(1 + (i%2)) * phiVals[i]**2
    C += (h/3)*phiVals[0] + (h/3)*phiVals[N]
    for i in range(len(phiVals)):
        phiVals[i] = phiVals[i]/np.sqrt(C)

def getPhiVals(epsilon, omega=1, xMax = 1, xMin=0, N=1000):
    phiVals = numerov(epsilon, omega=omega, xMax=xMax, xMin=xMin, N=N)
    normalizedPhiVals(phiVals, xMax=xMax, xMin=xMin, N=N)
    return phiVals

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
    # for i in range(len(lastPhiVals)):
    #     print(epsilonVals[i], lastPhiVals[i])
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

def shoot(epsilon1, epsilon2, omega, xMin, xMax, N, threshold=10**-3):
    epsilonHi = max(epsilon1, epsilon2)
    phiHi = numerov(epsilonHi, omega=omega, xMax = xMax, xMin=xMin, N=N)[-1]
    epsilonLo = min(epsilon1, epsilon2)
    phiLo = numerov(epsilonLo, omega=omega, xMax = xMax, xMin=xMin, N=N)[-1]
    epsilonMid = (epsilonHi + epsilonLo)/2
    phiMid = numerov(epsilonMid, omega=omega, xMax = xMax, xMin=xMin, N=N)[-1]
    while(abs(epsilonHi - epsilonLo) > threshold):
        # if phiLo * phiMid < 0:
        if (phiLo > 0 and phiMid < 0) or (phiLo < 0 and phiMid > 0):
            epsilonHi = epsilonMid
        else:
            epsilonLo = epsilonMid
        epsilonMid = (epsilonHi + epsilonLo)/2
        phiMid = numerov(epsilonMid, omega=omega, xMax = xMax, xMin=xMin, N=N)[-1]
        # print("%.4f %.4f" % (epsilonMid, phiMid))
    return epsilonMid

def getEigenvaluesByShooting(omega=1, xMax = 15, xMin=-15, N=1000, threshold=10**-3):
    numEigenvalues = 20 # remember to change this to 20
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
    prevPhiVal = numerov(epsilon, omega=omega, xMax = xMax, xMin=xMin, N=N)[-1]
    currPhiVal = None
    while(len(eigenvalues) < numEigenvalues):
        epsilon += epsilonStep
        currPhiVal = numerov(epsilon, omega=omega, xMax = xMax, xMin=xMin, N=N)[-1]
        # if we've encountered a sign change:
        # equivalent to currPhiVal * prevPhiVal < 0, but without overflow fears
        if ((currPhiVal > 0 and prevPhiVal < 0) or
                    (currPhiVal < 0 and prevPhiVal > 0)):
            shot = shoot(epsilon1=epsilon, epsilon2=epsilon - epsilonStep, 
                omega=omega, xMax = xMax, xMin=xMin, N=N, 
                threshold=threshold)
            eigenvalues.append(shot)
        prevPhiVal = currPhiVal
    return eigenvalues

def match(epsilon1, epsilon2, omega, xMin, xMax, N, threshold=10**-3):
    M = N//2




def foo():
    N = 10
    epsilonVals = [0.45 + (0.55-0.45)*i/N for i in range(N+1)]
    phiVals = [None for _ in range(len(epsilonVals))]
    for i in range(len(epsilonVals)):
        phiVals[i] = numerov(epsilonVals[i], omega=1, xMax = 11.0, xMin=-11.0, N=1000)[-1]
    try:
        import matplotlib.pyplot as plt 
        plt.plot(epsilonVals, phiVals)
        plt.show()
    except:
        pass
    for i in range(len(phiVals)):
        print(epsilonVals[i], phiVals[i])

def plotTrueWaveFunctions(eigenvalues, omega=1.0, xMin=-10, xMax=10, N=10000):
    #defining the following 2 because useful for plotting
    h = abs(xMax - xMin)/N
    xVals = [(xMin + i*h)*omega for i in range(N+1)]
    try:
        import matplotlib.pyplot as plt    
        for eigenvalue in eigenvalues:
            waveFunctionVals = getPhiVals(eigenvalue, omega=omega, xMax = xMax, xMin=xMin, N=N)
            plt.plot(xVals, waveFunctionVals)
        plt.show()
    except:
        pass


if __name__ == '__main__':
    # getWaveFunctionVals()
    # pass
    # foo()
    # e = shoot(epsilon1=0.4, epsilon2=0.55, 
    #     omega=1, xMin=-11, xMax=5, N=10000, threshold=10**-5)
    # print("e is", e)
    eigenvalues = getEigenvaluesByShooting(threshold=10**-10)
    print(eigenvalues)
    plotTrueWaveFunctions(eigenvalues)    
    # plotTrueWaveFunctions([(n + 0.5) for n in range(4)])
    
