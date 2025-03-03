## Lagrange interpolering
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import BarycentricInterpolator
        
def runge(x):
    return 1 / (1 + 25 * x**2)

def chebyshev(a, b, n):
    x = np.cos((2 * np.arange(1, n+1) - 1) * np.pi / (2 * n))
    return 0.5 * (a + b) + 0.5 * (b - a) * x

# Define the interval and number of points
a, b = -1, 1
n = 15

# Generate Chebyshev points
x_cheb = chebyshev(a, b, n)
y_cheb = runge(x_cheb)

# Use Chebyshev points for interpolation
inter = BarycentricInterpolator(x_cheb, y_cheb)

# Define a range of x values for plotting
x_plot = np.linspace(a, b, 1000)

# Calculate the Runge function and the interpolation
y_runge = runge(x_plot)
y_interpolated = inter(x_plot)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_plot, y_runge, label='Runge Function', color='blue')
plt.plot(x_plot, y_interpolated, label='Barycentric Interpolation', color='red', linestyle='--')
plt.scatter(x_cheb, y_cheb, color='black', zorder=5, label='Chebyshev Points')
plt.title('Runge Function and Barycentric Interpolation using Chebyshev Nodes')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()