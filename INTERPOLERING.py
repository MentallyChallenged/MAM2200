import matplotlib.pyplot as plt
import numpy as np

def funk(x):
    return 1 / (1 + 25 * x**2)

#x = np.linspace(-1, 1, 100)  
#plt.plot(x, y)  
#plt.show()  

a = -2
b = 3
n = 11

def chebyshev(a,b,n):
    x = np.cos((np.arange(1, n+1) * np.pi) / (n - 1))
    return 0.5*(a+b) + 0.5*(b-a)*x
print(chebyshev(a,b,n))

y = funk(chebyshev(a,b,n))  
#y = chebyshev(a,b,n)
plt.plot(chebyshev(a,b,n),y)
plt.show()
