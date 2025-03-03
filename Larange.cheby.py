## Lagrange interpolering
import numpy as np 
import matplotlib.pyplot as plt 

def vandermonde2(x, m):
    n = x.size
    V = np.ones((n,m))
    for i in range(1,m):
        V[:,i] = (x ** i) 
    return V
v = vandermonde2(np.array([[1,2,3]]), 4)

class juice:
    def __init__(self, x_verdier, y_verdier):
        self.n = x_verdier.size
        self.y = y_verdier
        VanderMatrise = vandermonde2(x_verdier, self.n)
        self.c = np.linalg.solve(VanderMatrise, y_verdier)
    def __call__(self, x):
        VanderMatrise = vandermonde2(x, self.n)
        return VanderMatrise @ self.c
        
def runge(x):
    return 1 / (1 + 25 * x**2)

def chebyshev(a, b, n):
    x = np.cos((2 * np.arange(1, n+1) - 1) * np.pi / (2 * n))
    return 0.5 * (a + b) + 0.5 * (b - a) * x

# Define the interval and number of points
a, b = -1, 1
n = 5

# Generate Chebyshev points
x_cheb = chebyshev(a, b, n)
y_cheb = runge(x_cheb)

# Use Chebyshev points for interpolation
inter = juice(x_cheb, y_cheb)

# Define a range of x values for plotting
x_plot = np.linspace(a, b, 1000)

# Calculate the Runge function and the interpolation
y_runge = runge(x_plot)
y_interpolated = inter(x_plot)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_plot, y_runge, label='Runge Function', color='blue')
plt.plot(x_plot, y_interpolated, label='Lagrange Interpolation', color='red', linestyle='--')
plt.scatter(x_cheb, y_cheb, color='black', zorder=5, label='Chebyshev Points')
plt.title('Runge Function and Lagrange Interpolation using Chebyshev Nodes')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()