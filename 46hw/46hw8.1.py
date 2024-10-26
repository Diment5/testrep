import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange, CubicSpline, PchipInterpolator, BPoly
def f(x):
    return 1 / (1 + x**2)
def f_prime(x):
    return -2 * x / (1 + x**2)**2
# Number of nodes
n_values = [5, 10, 15, 20]

# Interval
a, b = -5, 5

# Points to evaluate the interpolations
x_plot = np.linspace(a, b, 1000)
y_true = f(x_plot)
def lagrange_interpolation(n):
    x_nodes = np.linspace(a, b, n)
    y_nodes = f(x_nodes)
    # Construct Lagrange interpolating polynomial
    poly = lagrange(x_nodes, y_nodes)
    y_interp = poly(x_plot)
    return y_interp
def hermite_interpolation(n):
    x_nodes = np.linspace(a, b, n)
    y_nodes = f(x_nodes)
    dy_nodes = f_prime(x_nodes)
    # Prepare derivatives for Hermite interpolation
    derivatives = [(y, dy) for y, dy in zip(y_nodes, dy_nodes)]
    # Construct Hermite interpolating polynomial
    poly = BPoly.from_derivatives(x_nodes, derivatives)
    y_interp = poly(x_plot)
    return y_interp
def natural_cubic_spline(n):
    x_nodes = np.linspace(a, b, n)
    y_nodes = f(x_nodes)
    # Construct natural cubic spline
    cs = CubicSpline(x_nodes, y_nodes, bc_type='natural')
    y_interp = cs(x_plot)
    return y_interp
def clamped_cubic_spline(n):
    x_nodes = np.linspace(a, b, n)
    y_nodes = f(x_nodes)
    # Calculate derivatives at the endpoints
    fpa = f_prime(x_nodes[0])
    fpb = f_prime(x_nodes[-1])
    # Construct clamped cubic spline
    cs = CubicSpline(x_nodes, y_nodes, bc_type=((1, fpa), (1, fpb)))
    y_interp = cs(x_plot)
    return y_interp
methods = {
    'Lagrange': lagrange_interpolation,
    'Hermite': hermite_interpolation,
    'Natural Cubic Spline': natural_cubic_spline,
    'Clamped Cubic Spline': clamped_cubic_spline
}

for n in n_values:
    plt.figure(figsize=(14, 10))
    plt.plot(x_plot, y_true, 'k-', label='True Function', linewidth=2)
    for name, method in methods.items():
        try:
            y_interp = method(n)
            plt.plot(x_plot, y_interp, label=f'{name} (n={n})')
        except Exception as e:
            print(f"Error with {name} interpolation for n={n}: {e}")
    plt.title(f'Interpolation of f(x) with n={n} nodes')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.grid(True)
    plt.show()


