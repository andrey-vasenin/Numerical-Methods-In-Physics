import numpy as np

xmin = -8
xmax = 4
mesh = 1000

l = 1
n = 2
Z = 1

x, h = np.linspace(xmin, xmax, mesh)
r = np.exp(x) / Z
r1_2 = np.sqrt(r)
r2 = r ** 2
V = -1 / r

E = 0
Emin = V.min()
Emax = V.max()

ncount = -1

while ncount != n:
    E = (Emax + Emin) / 2

    f = 1 + (2 * r2 * (E - V) - (l + 1 / 2) ** 2) * h**2 / 12

    icl = V.searchsorted(E)

    # Solving inside the classically allowed region
    y[0] = 1
    y[1] = (6 - 5 * f[0]) * y[0] / f[1]
    for i in range(1, icl):
        y[i+1] = ((12 - 10 * f[i]) * y[i] - f[i-1] * y[i-1]) / f[i+1]
    yicl = y[icl]

    # Solving outside the classically forbidden region
    y[mesh] = h * (-1 if n % 2 > 0 else 1)
    y[mesh-1] = (12 - 10 * f[mesh]) * y[mesh] / f[mesh-1]
    for i in range(mesh-1, icl, -1):
        y[i-1] = ((12 - 10 * f[i]) * y[i] - f[i+1] * y[i+1]) / f[i-1]

    # Making the wave function continuous
        k = yicl / y[icl]
        y[icl:] *= k
    
    # Norming the solution
    y2 = y**2
    I = np.sum(y2[:-4] + 4*y2[2:-2]+ y2[4:]) * 2 * h / 3
    if I == 0:
        I = 0.1
    y = y / np.sqrt(I)

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
        if yLRder < 0:
            En = (Emax+temp)/2
            continue
        elif yLRder > 0:
            En = (Emin+temp)/2
            continue