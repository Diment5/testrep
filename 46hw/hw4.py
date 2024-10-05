from scipy.special import erf
import numpy as np
import matplotlib.pyplot as plt
from mypkg.Iteration1D import Iteration1D
alpha=0.138*10**-6
Ti =  20
Ts = -15
t=60**3 * 24
x = np.linspace(0,1,100)
y = erf(x/(2*(alpha*t)**1/2))*(Ti-Ts)+Ts 
plt.plot(x,y)
#plt.show()
f= lambda x:erf(x/(2*(alpha*t)**1/2))*(Ti-Ts)+Ts  
find3 = Iteration1D(f,'bisection')
find3.a = 0; find3.b = 1
find3.tol = 1e-13; find3.Nmax = 100
x3_bisection = find3.root()
print(x3_bisection)

find2 = Iteration1D(f,'newton')
find2.p0 = 1
find2.fp = lambda x: 23.35*np.exp(-0.35*x**2)
find2.tol = 1e-13; find2.Nmax = 100
x2_newton = find2.root()
print(x2_newton)

#secant 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

# Define the function and its derivative for Newton's method
def f(x):
    return x**6 - x - 1

def fp(x):
    return 6*x**5 - 1

# Use scipy to find the exact root of f(x)
alpha = root(f, 1.5).x[0]  # Find root using an initial guess of 1.5

# Secant method
def secant(f, x0, x1, tol, Nmax, alpha):
    e = []  # To store |x_k - alpha| values
    ier = 0
    count = 0
    
    if np.abs(f(x1) - f(x0)) < tol:
        ier = 1
        return e, count
    
    while count < Nmax:
        # Avoid division by zero
        if f(x1) == f(x0):
            ier = 1
            return e, count
        
        # Correct secant method formula
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        e.append(np.abs(x1 - alpha))  # Store |x_k - alpha|
        
        if np.abs(x2 - x1) < tol:
            break
        
        # Update x0 and x1 for the next iteration
        x0 = x1
        x1 = x2
        count += 1
    
    return e, count

# Newton method
def newton(f, fp, p0, tol, Nmax, alpha):
    e = []  # To store |x_k - alpha| values
    
    for it in range(Nmax):
        if fp(p0) == 0:
            print("Error: Derivative is zero.")
            return e, it
        
        p1 = p0 - f(p0) / fp(p0)
        e.append(np.abs(p0 - alpha))  # Store |x_k - alpha|
        
        if np.abs(p1 - p0) < tol:
            break
        
        p0 = p1
    
    return e, it

# Parameters
x0_secant = 1.0  # Initial guess 1 for secant method
x1_secant = 2.0  # Initial guess 2 for secant method
tol = 1e-6
Nmax = 20
p0_newton = 2.0  # Initial guess for Newton's method

# Run both methods
e_secant, count_secant = secant(f, x0_secant, x1_secant, tol, Nmax, alpha)
e_newton, count_newton = newton(f, fp, p0_newton, tol, Nmax, alpha)

# Compute the log-log plot data
e_secant_k1 = e_secant[1:]  # |x_{k+1} - alpha| for secant
e_secant_k = e_secant[:-1]  # |x_k - alpha| for secant

e_newton_k1 = e_newton[1:]  # |x_{k+1} - alpha| for newton
e_newton_k = e_newton[:-1]  # |x_k - alpha| for newton

# Plot the log-log graphs
plt.figure(figsize=(10,5))

# Secant method plot
plt.subplot(1, 2, 1)
plt.loglog(e_secant_k, e_secant_k1, 'o-', label='Secant Method')
plt.xlabel('$|x_k - \\alpha|$')
plt.ylabel('$|x_{k+1} - \\alpha|$')
plt.title('Secant Method: Log-Log Plot')
plt.legend()

# Newton method plot
plt.subplot(1, 2, 2)
plt.loglog(e_newton_k, e_newton_k1, 'o-', label='Newton Method')
plt.xlabel('$|x_k - \\alpha|$')
plt.ylabel('$|x_{k+1} - \\alpha|$')
plt.title('Newton Method: Log-Log Plot')
plt.legend()

plt.tight_layout()
plt.show()
