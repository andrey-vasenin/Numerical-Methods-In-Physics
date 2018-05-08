import numpy as np
import matplotlib.pyplot as plt

mesh = 10000 # number of intervals
xmin = 0
xmax = 10


x, h = np.linspace(xmin, xmax, num=mesh + 1, retstep=True) # grid
V = 150*(np.exp(-x+1) - 1)**2 # Potential
ones = np.ones(mesh + 1)

def solveMorsePot(n):

    Emin = V.min()
    Emax = V[mesh]
    y = np.zeros(mesh + 1) # solution
    ncount = -1 # Number of nodes

    yLRder = -1

    while Emax - Emin > 1e-6:
        En = (Emax + Emin) / 2
        f = 1 + (En - V) * h ** 2 / 6
        
        # Searching the right pivoting point
        imin = V.argmin()
        icl = V[imin:].searchsorted(En)+imin
        
        # Solving in the area before the right pivot point
        y[0] = 0
        y[1] = 1
        for i in range(1, icl):
            y[i+1] = ((12 - 10 * f[i]) * y[i] - f[i-1] * y[i-1]) / f[i+1]
        yicl = y[icl]

        # Solving in the right forbidden area (right -> left)
        y[mesh] = h if n % 2 > 0 else -h
        y[mesh - 1] = (12 - 10 * f[mesh]) * y[mesh] / f[mesh-1]
        for i in range(mesh-1, icl, -1):
            y[i-1] = ((12 - 10 * f[i]) * y[i] - f[i + 1] * y[i + 1]) / f[i - 1]
        
        # Making the solution continuous
        k = yicl / y[icl]
        y[icl:] *= k
        
        # Counting nodes
        ncount = 0
        sign = 1 if y[1] > 0 else -1
        for i in range(1, mesh):
            if y[i] < 0 and sign > 0:
                ncount += 1
                sign = -1
            elif y[i] > 0 and sign < 0:
                ncount += 1
                sign = 1
        
        # The first derivative
        yLRder = (y[icl + 1] + y[icl - 1] - (14 - 12 * f[icl]) * y[icl]) / h

        # Normalizing (Simpson's method)
        y2 = y**2
        I = np.sum(y2[:-4] + 4*y2[2:-2]+ y2[4:]) * 2 * h / 3
        if I == 0:
            I = 0.1
        y = y / np.sqrt(I)

        # Tuning the number of nodes
        if ncount < n:
            Emin = En
            continue
        elif ncount > n:
            Emax = En
            continue
        
        # Bisection method for sewing wave functions in classical and forbidden regions
        if n % 2 > 0:
            if yLRder < 0:
                Emax = En
                continue
            elif yLRder > 0:
                Emin = En
                continue
        else:
            if yLRder > 0:
                Emax = En
                continue
            elif yLRder < 0:
                Emin = En
                continue
        

    print("Energy of the " + str(n) + " level:", En)
    plt.fill_between(x, 5 * y + En, ones * En)

# Plotting the wave function
for N in range(0, 17):
    solveMorsePot(N)
plt.plot(x, V)
plt.grid(True)
plt.xlim(0, xmax)
plt.show()