import numpy as np
import matplotlib.pyplot as plt
# Skript som beregner interpolerende polynomer ved
# vandermonde-matrise og Lagrange-polynom :
# Definer funksjonen som vi skal integrere:
def f(x):
    return np.exp(-25*x**2)

# Funksjon som definerer Lagrange-polynomene:
def LagrangePol(k,xData,x):
    N=len(xData)
    L = np.ones(len(x))
    for j in range(N):
        if j != k:
            L = L*(x-xData[j])/(xData[k]-xData[j])
    return L

def chebyshev(a, b, n):
    x = np.cos((2 * np.arange(1, n+1) - 1) * np.pi / (2 * n))
    return 0.5 * (a + b) + 0.5 * (b - a) * x
a = -1
b = 1
N = 20
# Definer punkter som skal interpoleres:
xData = chebyshev(a,b,N)
# Chebyshev :
#xData = -np.cos(np.linspace(0,N-1,N)*np.pi/(N-1))
# y-koordinater til punktene:
yData = f(xData)

# Array for plotting av grafen til f:
x = np.linspace(-1,1,200)

# ---- Rett fram monomial interpolasjon ----
# Beregn vandermonde-matrisa og løs ligningssystemet:
va = np.vander(xData, len(xData))
pcoeffs = np.linalg.solve(va,yData)

# Definer polynomet på arrayet x:
pol = pcoeffs[0]*x**(N-1)
for k in range(N-1):
    pol += pcoeffs[k+1]*x**(N-1-(k+1))
    
    
# ---- Lagrange-polynom ----
# Definer polynomene L0, L1, ..., LN-1:
lpol = np.zeros(len(x))
for k in range(N):
    L = LagrangePol(k,xData,x)
    lpol += yData[k]*L
    
# Plott punkter og interpolerende polynom:
plt.plot(x,f(x))
plt.plot(xData,yData,'o','color','k')
plt.plot(x,lpol,'r')
plt.show()