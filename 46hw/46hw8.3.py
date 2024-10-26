import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.sin(10 * x)

def periodic_cubic_spline(x, y):
    n = len(x)
    h = np.diff(x)
    delta = np.diff(y) / h

    # Set up the cyclic tridiagonal system
    A = np.zeros((n, n))
    rhs = np.zeros(n)

    # Coefficients for the tridiagonal part
    for i in range(1, n - 1):
        A[i, i - 1] = h[i - 1]
        A[i, i] = 2 * (h[i - 1] + h[i])
        A[i, i + 1] = h[i]
        rhs[i] = 6 * (delta[i] - delta[i - 1])

    # Adjust the first and last equations for periodicity
    # First equation
    A[0, 0] = 2 * (h[0] + h[-1])
    A[0, 1] = h[0]
    A[0, -1] = h[-1]
    rhs[0] = 6 * (delta[0] - delta[-1])

    # Last equation
    A[-1, -2] = h[-2]
    A[-1, -1] = 2 * (h[-2] + h[-1])
    A[-1, 0] = h[-1]
    rhs[-1] = 6 * (delta[-1] - delta[-2])

    # Solve the system
    m = np.linalg.solve(A, rhs)

    # Compute spline coefficients
    spline_coeffs = []
    for i in range(n - 1):
        hi = h[i]
        mi = m[i]
        mi1 = m[i + 1]
        yi = y[i]
        yi1 = y[i + 1]

        a = (mi1 - mi) / (6 * hi)
        b = mi / 2
        c = (yi1 - yi) / hi - (2 * hi * mi + hi * mi1) / 6
        d = yi

        spline_coeffs.append((a, b, c, d))

    return spline_coeffs

def evaluate_spline(x, spline_coeffs, xi):
    # Find the interval xi is in
    n = len(spline_coeffs)
    if xi < x[0]:
        i = 0
    elif xi > x[-1]:
        i = n - 1
    else:
        i = np.searchsorted(x, xi) - 1

    a, b, c, d = spline_coeffs[i]
    hi = xi - x[i]
    si = a * hi**3 + b * hi**2 + c * hi + d
    return si

# Number of nodes
n = 20
x_nodes = np.linspace(0, 2 * np.pi, n + 1)[:-1]  # Exclude the last point to make it periodic
y_nodes = f(x_nodes)

# Compute the periodic cubic spline coefficients
spline_coeffs = periodic_cubic_spline(x_nodes, y_nodes)

# Evaluate the spline at a set of points
x_eval = np.linspace(0, 2 * np.pi, 1000)
y_eval = np.array([evaluate_spline(x_nodes, spline_coeffs, xi % (2 * np.pi)) for xi in x_eval])

# Plot the function and the spline approximation
plt.figure(figsize=(10, 6))
plt.plot(x_eval, f(x_eval), label='True Function')
plt.plot(x_eval, y_eval, label='Periodic Cubic Spline', linestyle='--')
plt.scatter(x_nodes, y_nodes, color='red', label='Nodes')
plt.legend()
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Periodic Cubic Spline Interpolation of sin(10x)')
plt.grid(True)
plt.show()
