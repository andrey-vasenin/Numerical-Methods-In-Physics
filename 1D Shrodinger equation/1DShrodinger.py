import numpy as np
import matplotlib.pyplot as plt

mesh = 2000
xmax = 10
n = 0
n2 = (n - n % 2)/2
En = 1/2
yLRder_E = 1

x, h = np.linspace(0, xmax, num=mesh+1, retstep=True)
y = np.empty(mesh + 1)
f = np.empty(mesh + 1)

En = n + 0.5
f[0] = 1 + En * h**2 / 6
f[1] = 1 + (2 * En - x[1]**2) * h**2 / 12
if n % 2 > 0:
    y[0] = 0
    y[1] = 1
else:
    y[0] = 1
    y[1] = (6 - 5 * f[0]) * y[0] / f[1]

icl = mesh
for i in range(0, mesh):
    if En < x[i]**2/2:
        icl = i
        break
y[icl] = 0

for i in range(1, icl):
    f[i+1] = 1 + (2 * En - x[i+1]**2) * h**2 / 12
    y[i+1] = ((12 - 10 * f[i]) * y[i] - f[i-1] * y[i-1]) / f[i+1]
yicl = y[icl]

# Outside the classical region
y[mesh] = h
f[mesh] = 1 + (2 * En - x[mesh]**2) * h**2 / 12
f[mesh-1] = 1 + (2 * En - x[mesh-1]**2) * h**2 / 12
y[mesh-1] = (12 - 10 * f[mesh]) * y[mesh] / f[mesh-1]
for i in range(mesh-1, icl, -1):
    f[i-1] = 1 + (2 * En - x[i-1]**2) * h**2 / 12
    y[i-1] = ((12 - 10 * f[i]) * y[i] - f[i+1] * y[i+1]) / f[i-1]

# Making the wave function continuous
k = yicl / y[icl]
for i in range(icl, mesh+1):
    y[i] *= k

# Comparing first derivatives
yLRder = (y[icl+1] + y[icl-1] - (14 - 12 * f[icl]) * y[icl])/h

# Normalizing
I = 0
i = 0
for i in range(0, mesh-4):
    I += (y[i] + 4*y[i+2] + y[i+4])*2*h/3
print(I)
for i in range(0, mesh+1):
    y[i] /= I
I = 0
i = 0
for i in range(0, mesh-4):
    I += (y[i] + 4*y[i+2] + y[i+4])*2*h/3
print(I)
print(En, x[icl])
plt.plot(x, y, "-b")
plt.plot(-x, (-1)**n * y, "-b")
plt.plot([x[icl], x[icl]], [-0.4, 0.4], ":r")
plt.grid(True)
plt.show()