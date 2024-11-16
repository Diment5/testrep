import numpy as np

def integral_trapezoidal(f, a, b, n):
    """
    Composite Trapezoidal Rule
    :param f: Function to integrate
    :param a: Start of the interval
    :param b: End of the interval
    :param n: Number of subintervals
    :return: Approximation of the integral
    """
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = f(x)
    integral = h * (0.5 * y[0] + np.sum(y[1:-1]) + 0.5 * y[-1])
    return integral

def integral_simpson(f, a, b, n):
    """
    Composite Simpson's Rule
    :param f: Function to integrate
    :param a: Start of the interval
    :param b: End of the interval
    :param n: Number of subintervals (must be even)
    :return: Approximation of the integral
    """
    if n % 2 != 0:
        raise ValueError("Number of subintervals (n) must be even for Simpson's rule.")
    
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = f(x)
    integral = h / 3 * (y[0] + 4 * np.sum(y[1:n:2]) + 2 * np.sum(y[2:n-1:2]) + y[-1])
    return integral

def compute_integral(method="trapezoidal", n=10):
    """
    Compute the integral using the specified method
    :param method: 'trapezoidal' or 'simpson'
    :param n: Number of subintervals
    :return: Approximation of the integral
    """
    # Define the function
    f = lambda s: 1 / (1 + s**2)
    a, b = -5, 5  # Interval bounds

    if method == "trapezoidal":
        return integral_trapezoidal(f, a, b, n)
    elif method == "simpson":
        return integral_simpson(f, a, b, n)
    else:
        raise ValueError("Unknown method. Choose 'trapezoidal' or 'simpson'.")

# Example usage
n_trapezoidal = 10  # Number of subintervals for Trapezoidal rule
n_simpson = 10      # Number of subintervals for Simpson's rule (must be even)

result_trapezoidal = compute_integral(method="trapezoidal", n=n_trapezoidal)
result_simpson = compute_integral(method="simpson", n=n_simpson)

print(f"Trapezoidal Rule Result (n={n_trapezoidal}): {result_trapezoidal:.6f}")
print(f"Simpson's Rule Result (n={n_simpson}): {result_simpson:.6f}")
