import numpy as np
import matplotlib.pyplot as plt

mesh = 10000 # number of intervals
xmax = 10

x, h = np.linspace(0, xmax, num=mesh+1, retstep=True) # grid
V = x**2/2 # Potential
ones = np.ones(2 * mesh + 2)

def solveQuantOsc(n):
    n2 = (n - n % 2)/2 # Number of nodes on the rigth side

    y = np.zeros(mesh + 1) # solution

    Emin = V.min()
    Emax = V.max()

    ncount = -1 # Number of nodes

    yLRder = -1 # y'(left) - y'(right)
    yLRder_prev = -1
    En = (Emax + Emin)/2 # Energy of the n-th level
    En_prev = 0

    while ncount != n2 or np.abs(yLRder) > 1e-7:
        
        f = 1 + (En - V) * h**2 / 6
        icl = V.searchsorted(En) # The index of the beginning of classically fobidden area

        # Solving inside the classical region (left -> right)
        y[0] = 0 if n % 2 > 0 else 1
        y[1] = 1 if n % 2 > 0 else (6 - 5 * f[0]) * y[0] / f[1]
        for i in range(1, icl):
            y[i+1] = ((12 - 10 * f[i]) * y[i] - f[i-1] * y[i-1]) / f[i+1]
        yicl = y[icl]

        # Solving outside the classical region (right -> left)
        y[mesh] = h
        y[mesh-1] = (12 - 10 * f[mesh]) * y[mesh] / f[mesh-1]
        for i in range(mesh-1, icl, -1):
            y[i-1] = ((12 - 10 * f[i]) * y[i] - f[i+1] * y[i+1]) / f[i-1]

        # Making the wave function continuous
        k = yicl / y[icl]
        y[icl:] *= k

        # Normalizing (Simpson's method)
        y2 = y**2
        I = np.sum(y2[:-4] + 4*y2[2:-2]+ y2[4:]) * 2 * h / 3
        if I == 0:
            I = 0.1
        # print(np.sqrt(I))
        y = y / np.sqrt(I)
        
        # Uncomment to see intermediate plots
        # plt.plot(x, y, "-b")
        # plt.plot(-x, (-1)**n * y, "-b")
        # plt.plot([x[icl], x[icl]], [-0.4, 0.4], ":r")
        # plt.grid(True)
        # plt.show()
        
        # Counting nodes
        ncount = 0
        sign = 1 if y[1] > 0 else -1
        for i in range(1, icl+1):
            if y[i] < 0 and sign > 0:
                ncount += 1
                sign = -1
            elif y[i] > 0 and sign < 0:
                ncount += 1
                sign = 1
        
        # First derivative
        yLRder_prev = yLRder
        yLRder = (y[icl+1] + y[icl-1] - (14 - 12 * f[icl]) * y[icl])/h

        # Tuning the energy
        if ncount < n2:
            En_prev = En
            Emin = En
            En = (Emax+Emin)/2
            continue
        elif ncount > n2:
            En_prev = En
            Emax = En
            En = (Emax+Emin)/2
            continue
        
        # # Newton's method for sewing wave functions in classical and forbidden regions
        temp = En
        En = En - yLRder / (yLRder - yLRder_prev) * (En - En_prev)
        En_prev = temp
        # Use bisection method in case Newton's method behaves badly
        if En > Emax or En < Emin:
            # if n % 2 > 0:
            if yLRder < 0:
                # Emax = temp
                # En = (Emax+Emin)/2
                En = (Emax+temp)/2
                continue
            elif yLRder > 0:
                # Emin = temp
                # En = (Emax+Emin)/2
                En = (Emin+temp)/2
                continue
            # else:
            #     if yLRder > 0:
            #         Emax = temp
            #         En = (Emax+Emin)/2
            #         continue
            #     elif yLRder < 0:
            #         Emin = temp
            #         En = (Emax+Emin)/2
            #         continue

    print("Energy of the " + str(n) + " level:", En)
    plt.fill_between(np.append(-x[::-1], x), np.append((-1)**n*y[::-1], y)+En, ones * En)

# Plotting the wave function
# plt.plot(x, y, "-b")
# plt.plot(-x, (-1)**n * y, "-b")
# # plt.plot([x[icl], x[icl]], [-0.4, 0.4], ":r")
for N in range(0, 10):
    solveQuantOsc(N)
plt.plot(np.append(-x[::-1], x), np.append(V[::-1], V))
plt.grid(True)
plt.xlim(-xmax, xmax)
plt.show()