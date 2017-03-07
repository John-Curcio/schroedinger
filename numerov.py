"""
LANG: Python 2
Use Numerov method, together with either shooting method or matching method, for 
solving energy eigenvalues of 1D Schroedinger equation for electron (some eqn)
"""
omega = 1
N = 20
h = 0.05
x_min = 0
x_max = x_min * h*N

def get_phi(phi_hat, h):
    phi_hat = phi/(1 - (h**2)/12)

def get_phi_hat(phi, h):
    phi_hat = (1 - (h**2)/12)*phi

def f(x, epsilon):
    return 0.25 * (omega**2 ) * x**2 - epsilon

def get_phi_2(x_1, phi_1, phi_0, epsilon):
    phi_0_hat = get_phi_hat(phi_0, h)
    phi_1_hat = get_phi_hat(phi_1, h)
    phi_2_hat = 2*phi_1_hat + (h**2)*f(x_1, epsilon)*phi_1 - phi_0_hat
    return get_phi(phi_2_hat)

def get_x_vals():
    x_vals = [None]*N
    for i in xrange(0, N):
        x_vals[i] = i*h + x_min
    return x_vals

def shoot():
    x_vals = get_x_vals()
    eigenstates = [1]*N
    eigenvalues = [1]*N
    eigenstates[0] = 0
    eigenstates[1] = 1
    for i in xrange(2, N):
        while(not almost_eq(eigenstates[i]**2, 0))
            eigenstates[i] = get_phi_2(x_vals[i-1], eigenstates[i-1], eigenstates[i-2], eigenvalues[i-1])

if __name__ == '__main__':
    print f(1, 1)
    #main()