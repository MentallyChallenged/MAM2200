import numpy as np
import matplotlib.pyplot as plt

def f(xs):
#return xs**2 + 0.2*np.sin(10*xs)
    return np.exp(xs-np.sqrt(xs))-xs

def fd(xs):
#return 2*xs + 0.2*10*np.cos(10*xs)
    return np.exp(xs-np.sqrt(xs))*(1-1/(2*np.sqrt(xs)))-1

def descent(f,fd,xp,gamma=1,kmax=100,tol=1.e-8):
    for k in range(1,kmax):
        xn = xp - gamma*fd(xp)
        err = np.abs(xn-xp)
        print(k, xn, err, f(xn))
        plt.plot([xp,xn],[f(xp),f(xn)],color='red')

        if err < tol:
            break
        xp = np.copy(xn)
    return xn
# Startverdi:
xp = 3
xn = descent(f,fd,xp)
print(xn)
xx=np.linspace(0,3,200)
plt.plot(xx,f(xx),color='black')
plt.show()
