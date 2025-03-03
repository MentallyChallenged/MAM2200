## Kubisk spline interpolasjon 

import numpy as np
import matplotlib.pyplot as plt
# Definer funksjonen som vi skal integrere:
def f(x):
#return x**2*np.exp(-x)
    return np.exp(-np.abs(x))

def computecs(dataxs,datays):
    n = dataxs.size
    A = np.zeros((n-2,n-2))
    np.fill_diagonal(A, 2*(dataxs[2:]-dataxs[:-2]))
    np.fill_diagonal(A[1:,:], dataxs[2:-1]-dataxs[1:-2])
    np.fill_diagonal(A[:,1:], dataxs[2:-1]-dataxs[1:-2])
    b1 = (datays[2:]-datays[1:-1])/(dataxs[2:]-dataxs[1:-1])
    b2 = (datays[1:-1]-datays[:-2])/(dataxs[1:-1]-dataxs[:-2])
    bs = 6*(b1 - b2)
    cs = np.zeros(n)
    cs[1:-1] = np.linalg.solve(A,bs)
    return cs
def splineinterp(dataxs,datays,cs,x):
    k = np.argmax(dataxs>x)
    xk = dataxs[k]; xk1 = dataxs[k-1]
    yk = datays[k]; yk1 = datays[k-1]
    ck = cs[k]; ck1 = cs[k-1]
    val = yk1*(xk-x)/(xk-xk1) + yk*(x-xk1)/(xk-xk1)
    val -= ck1*((xk-x)*(xk-xk1) - (xk-x)**3/(xk-xk1))/6
    val -= ck*((x-xk1)*(xk-xk1) - (x-xk1)**3/(xk-xk1))/6
    return val

# Data som skal interpoleres
N = 16
xData = np.linspace(-1,1,N)
yData = f(xData)

# Regn ut den kubiske splinen:
cs = computecs(xData, yData)
M = 256
x=np.linspace(-1,1,M)

#x = 0.95;
pofx = np.zeros(M)
for k in range(M):
    pofx[k] = splineinterp(xData, yData, cs, x[k])

plt.plot(xData, yData,'+r',markersize=12)
plt.plot(x,pofx)
plt.show()
