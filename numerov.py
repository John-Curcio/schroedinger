import numpy as np
import sys

# LANG: Python 3
# Use Numerov method, together with either shooting method or matching method, for 
# solving energy eigenvalues of 1D Schroedinger equation for electron (some eqn)

def externalPotential(x, omega=1):
    return 0.25 * (omega**2) * (x**2)

def f(epsilon, x, omega=1):
    return externalPotential(x, omega=omega) - epsilon

def numerov(epsilon, omega=1, xMax = 1, xMin=0, N=1000):
    h = (xMax - xMin)/N
    xVals = [xMin + i*h for i in range(N+1)]
    assert(abs(xVals[-1] - xMax) < 10**-5)
    fVals = [f(epsilon, xVals[i], omega) for i in range(N+1)]
    phiVals = [None]*(N+1)
    phiVals[0] = 0
    phiVals[1] = h
    for i in range(1, N):
        prevPhiHat = (1 - (1/12)*(h**2)*fVals[i-1])*phiVals[i-1]
        nextPhiHat = (2 + (5/6)*(h**2)*fVals[i])*phiVals[i] - prevPhiHat
        phiVals[i+1] = nextPhiHat / (1 - (1/12)*(h**2)*fVals[i+1])
    return phiVals

def normalizedPhiVals(phiVals, xMax, xMin, N):
    h = (xMax - xMin)/N
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
    h = (xMax - xMin)/N
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
        # equivalent to "if phiLo * phiMid < 0:", but without fear of overflow
        if (phiLo > 0 and phiMid < 0) or (phiLo < 0 and phiMid > 0):
            epsilonHi = epsilonMid
            phiHi = phiMid
        else:
            epsilonLo = epsilonMid
            phiLo = phiMid
        epsilonMid = (epsilonHi + epsilonLo)/2
        phiMid = numerov(epsilonMid, omega=omega, xMax = xMax, xMin=xMin, N=N)[-1]
        # print("%.4f %.4f" % (epsilonMid, phiMid))
    return epsilonMid

def match(epsilon1, epsilon2, omega, xMin, xMax, N, threshold=10**-3):
    epsilonHi = max(epsilon1, epsilon2)
    epsilonLo = min(epsilon1, epsilon2)
    epsilonMid = (epsilonLo + epsilonHi)/2
    if abs(epsilonHi - epsilonLo) < threshold:
        return epsilonMid
    deltaHi = getDelta(epsilonHi, omega=omega, xMin=xMin, xMax=xMax, N=N)
    deltaMid = getDelta(epsilonMid, omega=omega, xMin=xMin, xMax=xMax, N=N)
    deltaLo = getDelta(epsilonLo, omega=omega, xMin=xMin, xMax=xMax, N=N)
    maxIters = 500
    iters = 0
    while(abs(epsilonHi - epsilonLo) > threshold and iters < maxIters):
        iters += 1
        if (deltaHi < 0 and deltaMid > 0) or (deltaHi > 0 and deltaMid < 0):
            deltaLo = deltaMid
            epsilonLo = epsilonMid
        else:
            deltaHi = deltaMid #not necessary, but w/e
            epsilonHi = epsilonMid
        epsilonMid = (epsilonLo + epsilonHi)/2
        deltaMid = getDelta(epsilonMid, omega=omega, xMin=xMin, xMax=xMax, N=N)
    return epsilonMid

def getEigenvalues(method=shoot, omega=1, xMax = 15, xMin=-15, N=1000, threshold=10**-3, numEigenvalues=20):
    epsilon = 0 
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
            eigenvalue = method(epsilon1=epsilon, epsilon2=epsilon - epsilonStep, 
                omega=omega, xMax = xMax, xMin=xMin, N=N, 
                threshold=threshold)
            eigenvalues.append(eigenvalue)
        prevPhiVal = currPhiVal
    return eigenvalues


def getDelta(epsilon=1.5, omega=1, xMin=-10.0, xMax=10.0, N=10000):
    M = N//2
    h = (xMax - xMin)/N
    phiLeftVals = numerov(epsilon, omega=omega, xMax = xMax/2, xMin=xMin, N=N//2)
    phiRightVals = numerov(epsilon, omega=omega, xMax = xMin, xMin=xMax/2, N=N//2)
    # scaling - don't actually need to normalize the other values, since
    # we don't use them to compute delta
    CLeft = (phiLeftVals[-1])
    phiLeftVals[-1] /= CLeft
    phiLeftVals[-2] /= CLeft
    CRight = (phiRightVals[-1])
    phiRightVals[-1] /= CRight
    phiRightVals[-2] /= CRight
    delta = phiLeftVals[-2] + phiRightVals[-2] 
    delta -= (2 + f(epsilon, xMax/2, omega=omega)*h**2)*phiLeftVals[-1]
    delta /= h
    return delta

def plotTrueWaveFunctions(eigenvalues, omega=1.0, xMin=-10, xMax=10, N=10000, title="Plot of Numerov Calculation for Given Epsilon"):
    #defining the following 2 because useful for plotting
    h = (xMax - xMin)/N
    xVals = [(xMin + i*h)*omega for i in range(N+1)]
    try:
        import matplotlib.pyplot as plt    
        for eigenvalue in eigenvalues:
            waveFunctionVals = getPhiVals(eigenvalue, omega=omega, xMax = xMax, xMin=xMin, N=N)
            plt.plot(xVals, waveFunctionVals)
        plt.ylim(-1, 1)
        plt.ylabel("Amplitude of Wave Function $\phi$")
        plt.title(title)
        plt.xlim(xMin, xMax)
        plt.xlabel("Positional Parameter $x$")
        plt.show()
    except:
        print("couldn't plot for some reason")


if __name__ == '__main__':
    omega = float(input("\nspring coefficient (omega): "))
    N = int(input("number of uniform mesh points (N): "))
    xMin = float(input("left boundary condition (xMin): "))
    xMax = float(input("right boundary condition (xMax): "))
    threshold = float(input("precision threshold for eigenvalue calculation (e.g. 0.0001): "))
    numEigenvalues = int(input("how many of the first energy eigenvalues to calculate? (>=1): "))
    plot = input("Would you like to plot the resulting wavefunctions? (Y/N): ")
    plot = plot.strip()
    numericalEigenvalues = getEigenvalues(method=shoot, omega=omega, 
        xMax = xMax, xMin=xMin, N=N, threshold=threshold, numEigenvalues=numEigenvalues)
    analyticalEigenvalues = [(n + 0.5)*omega for n in range(numEigenvalues)]
    print("First 20 Eigenvalues:\nAnalytical | Numerical")
    for i in range(numEigenvalues):
        print("%.4f     | %.4f"%(analyticalEigenvalues[i], numericalEigenvalues[i]))
    if plot == "Y" or plot == "y":
        plotTrueWaveFunctions(numericalEigenvalues, omega=omega, xMin=xMin, 
            xMax=xMax, N=N, 
            title="Plot of Wavefunction with Numerically-Calculated Eigenvalues")
        plotTrueWaveFunctions(analyticalEigenvalues, omega=omega, xMin=xMin, 
            xMax=xMax, N=N, 
            title="Plot of Wavefunction with Analytically-Calculated Eigenvalues")
