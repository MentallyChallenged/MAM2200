## Python er klasse orientert btw 
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

x = np.linspace(-1, 1, 5)
y = runge(x)

inter = juice(x,y)
# Define a range of x values for plotting
x_plot = np.linspace(-1, 1, 1000)

# Calculate the Runge function and the interpolation
y_runge = runge(x_plot)
y_interpolated = inter(x_plot)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_plot, y_runge, label='Runge Function', color='blue')
plt.plot(x_plot, y_interpolated, label='Lagrange Interpolation', color='red', linestyle='--')
plt.scatter(x, y, color='black', zorder=5, label='Interpolation Points')
plt.title('Runge Function and Lagrange Interpolation')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()