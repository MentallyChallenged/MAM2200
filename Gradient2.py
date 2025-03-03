## Gradient decent 

import numpy as np

#Her gir vi orginal funksjonen vår
def funksjon(Xs):
    x0, x1 = Xs
    return x0**2 - 2*x0 + x1**4 - 2*x1**2 + x1 

#Her gir gradienten av funksjonen vår, Så når vi har partiel derviert for per variabel 
def funskjon_gradientet(Xs):
    x0,x1 = Xs 
    return np.array([2*x0 -2 , 4*x1**3 - 4*x1 + 1])

def sum_av_differensen(x_old,x_new):
    error = np.sum(np.abs((x_new - x_old)/x_new))
    return error


def desent(
    phi,gradient, x_old,gamma = 0.15, #steglengde,
    K_max = 200, # maks antall iterasjon
    tolaranse = 0.00000001):
    
    for k in range(1,K_max):
        x_new = x_old - gamma*gradient(x_old)
        error = sum_av_differensen(x_old, x_new)
        print(k, x_new, phi(x_new))
        if error <= tolaranse:
            break
        x_old = np.copy(x_new)
    return x_new

bunnpunkt = desent(funksjon, funskjon_gradientet, np.array([2,0.25]))
print("BUNNPUNKT: ", bunnpunkt)

a = np.array([1.0,2.0,3.0,4.0])
n = a.size
h = 0.00001
a_pluss = a

for i in range(n):
    a_pluss = a.copy()
    a_pluss[i] += h
    print(i, a_pluss)