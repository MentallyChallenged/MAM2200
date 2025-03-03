import numpy as np

def vandermonde(x):
    n = x.size
    V = np.zeros((n, n))
    for i in range(n):
        V[:, i] = x ** i
    return V
m = vandermonde(np.array([[1, 2, 3]]))
print(m)

def vandermonde2(x, m):
    n = x.size
    V = np.ones((n,m))
    for i in range(1,m):
        V[:,i] = (x ** i) 
    return V
d = vandermonde2(np.array([[0,1,2,3]]), 4)
print(d)