import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Sub-question (a): Both Numerator and Denominator are Cubic Polynomials [3/3]
def pade_approximation_a():
    print("Sub-question (a): Both Numerator and Denominator are Cubic Polynomials [3/3]")
    # Define the symbolic variable
    x = sp.symbols('x')

    # Define symbolic coefficients for P(x) and Q(x)
    p0, p1, p2, p3 = sp.symbols('p0 p1 p2 p3')
    q1, q2, q3 = sp.symbols('q1 q2 q3')

    # Define P(x) and Q(x)
    P = p0 + p1*x + p2*x**2 + p3*x**3
    Q = 1 + q1*x + q2*x**2 + q3*x**3

    # Define the Taylor series expansion of sin(x) up to x^7
    sinx_series = x - x**3/6 + x**5/120

    # Multiply Q(x) and sin(x), expand, and collect terms up to x^6
    Q_sin = Q * sinx_series
    Q_sin = sp.series(Q_sin, x, 0, 7).removeO()
    Q_sin = sp.collect(Q_sin, x)

    # Compute P(x) - Q(x) * sin(x)
    expr = P - Q_sin

    # Collect coefficients of like powers of x
    coeffs = sp.collect(expr, x, evaluate=False)

    # Set up equations by equating coefficients to zero
    equations = []
    for power in range(7):  # From x^0 to x^6
        coeff = coeffs.get(x**power, 0)
        equations.append(sp.Eq(coeff, 0))

    # Solve the system of equations
    variables = [p0, p1, p2, p3, q1, q2, q3]
    solutions = sp.solve(equations, variables)

    # Display the solutions
    print("Solutions for coefficients:")
    for var in variables:
        print(f"{var} = {solutions.get(var)}")
    print("\n")

# Sub-question (b): Numerator is Quadratic and Denominator is Quartic [2/4]
def pade_approximation_b():
    print("Sub-question (b): Numerator is Quadratic and Denominator is Quartic [2/4]")
    # Define the symbolic variable
    x = sp.symbols('x')

    # Define symbolic coefficients for P(x) and Q(x)
    p0, p1, p2 = sp.symbols('p0 p1 p2')
    q1, q2, q3, q4 = sp.symbols('q1 q2 q3 q4')

    # Define P(x) and Q(x)
    P = p0 + p1*x + p2*x**2
    Q = 1 + q1*x + q2*x**2 + q3*x**3 + q4*x**4

    # Define the Taylor series expansion of sin(x) up to x^7
    sinx_series = x - x**3/6 + x**5/120

    # Multiply Q(x) and sin(x), expand, and collect terms up to x^6
    Q_sin = Q * sinx_series
    Q_sin = sp.series(Q_sin, x, 0, 7).removeO()
    Q_sin = sp.collect(Q_sin, x)

    # Compute P(x) - Q(x) * sin(x)
    expr = P - Q_sin

    # Collect coefficients of like powers of x
    coeffs = sp.collect(expr, x, evaluate=False)

    # Set up equations by equating coefficients to zero
    equations = []
    for power in range(7):  # From x^0 to x^6
        coeff = coeffs.get(x**power, 0)
        equations.append(sp.Eq(coeff, 0))

    # Solve the system of equations
    variables = [p0, p1, p2, q1, q2, q3, q4]
    solutions = sp.solve(equations, variables)

    # Display the solutions
    print("Solutions for coefficients:")
    for var in variables:
        print(f"{var} = {solutions.get(var)}")
    print("\n")

# Sub-question (c): Numerator is Quartic and Denominator is Quadratic [4/2]
def pade_approximation_c():
    print("Sub-question (c): Numerator is Quartic and Denominator is Quadratic [4/2]")
    # Define the symbolic variable
    x = sp.symbols('x')

    # Define symbolic coefficients for P(x) and Q(x)
    p0, p1, p2, p3, p4 = sp.symbols('p0 p1 p2 p3 p4')
    q1, q2 = sp.symbols('q1 q2')

    # Define P(x) and Q(x)
    P = p0 + p1*x + p2*x**2 + p3*x**3 + p4*x**4
    Q = 1 + q1*x + q2*x**2

    # Define the Taylor series expansion of sin(x) up to x^7
    sinx_series = x - x**3/6 + x**5/120

    # Multiply Q(x) and sin(x), expand, and collect terms up to x^6
    Q_sin = Q * sinx_series
    Q_sin = sp.series(Q_sin, x, 0, 7).removeO()
    Q_sin = sp.collect(Q_sin, x)

    # Compute P(x) - Q(x) * sin(x)
    expr = P - Q_sin

    # Collect coefficients of like powers of x
    coeffs = sp.collect(expr, x, evaluate=False)

    # Set up equations by equating coefficients to zero
    equations = []
    for power in range(7):  # From x^0 to x^6
        coeff = coeffs.get(x**power, 0)
        equations.append(sp.Eq(coeff, 0))

    # Solve the system of equations
    variables = [p0, p1, p2, p3, p4, q1, q2]
    solutions = sp.solve(equations, variables)

    # Display the solutions
    print("Solutions for coefficients:")
    for var in variables:
        print(f"{var} = {solutions.get(var)}")
    print("\n")

# Main execution
if __name__ == "__main__":
    pade_approximation_a()
    pade_approximation_b()
    pade_approximation_c()

# Define the interval
x = np.linspace(0, 5, 500)

# True function
sin_x = np.sin(x)

# Padé Approximant (a) and (c): [3/3] and [4/2]
def pade_ac(x):
    numerator = x - (7/60)*x**3
    denominator = 1 + (1/20)*x**2
    return numerator / denominator

# Padé Approximant (b): [2/4]
def pade_b(x):
    denominator = 1 + (1/6)*x**2 + (7/360)*x**4
    return x / denominator

# Sixth-order Maclaurin Polynomial
def maclaurin(x):
    return x - (x**3)/6 + (x**5)/120

# Calculate approximations
pade_ac_values = pade_ac(x)
pade_b_values = pade_b(x)
maclaurin_values = maclaurin(x)

# Calculate errors
error_pade_ac = np.abs(sin_x - pade_ac_values)
error_pade_b = np.abs(sin_x - pade_b_values)
error_maclaurin = np.abs(sin_x - maclaurin_values)

# Plotting the errors
plt.figure(figsize=(12, 6))
plt.plot(x, error_pade_ac, label='Padé Approximation a) and c)')
plt.plot(x, error_pade_b, label='Padé Approximation b)')
plt.plot(x, error_maclaurin, label='Maclaurin')
plt.title('Error Comparison')
plt.xlabel('x')
plt.ylabel('Absolute Error')
plt.legend()
plt.grid(True)
plt.show()