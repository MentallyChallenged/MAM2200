import matplotlib.pyplot as plt
import numpy as np

def funk(x):
    return np.sqrt(4 - x**2)

a = -2
b = 2
n = 20

def chebyshev(a, b, n):
    x = np.cos((np.arange(1, n+1) * np.pi) / (n - 1))
    return 0.5 * (a + b) + 0.5 * (b - a) * x

# Generer Chebyshev-punktene én gang
x_cheb = chebyshev(a, b, n)
y_cheb = funk(x_cheb)

# Plotter punktene
plt.plot(x_cheb, y_cheb, "o", label="Chebyshev-punkter")
plt.plot(x_cheb, y_cheb, "x", label="F(x) verdier")

# Tegner funksjonen for sammenligning
x_vals = np.linspace(a, b, 100)
y_vals = funk(x_vals)
plt.plot(x_vals, y_vals, label="f(x)")

plt.axis("equal")
plt.grid()
plt.legend()
plt.show()
