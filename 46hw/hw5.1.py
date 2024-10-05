import numpy as np
import math
from numpy.linalg import norm 

x0=[1,1]
def evalF(x): 

    F = np.zeros(2)
    
    F[0] = 3*x[0]**2-x[1]**2
    F[1] = 3*x[0]*x[1]**2-x[0]**3-1
    
    return F

Jinv1 = np.array([0]*4).reshape(2,2)
Jinv1[0,:]=[1/6,1/18]
Jinv1[1,:]=[0,1/6]


tol1=10**-14

Nmax1=30
count=0
def iterate(x0,Nmax,Jinv):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    for its in range(Nmax):
       
       F = evalF(x0)
       
       x1 = x0 - Jinv1.dot(F)
       
       
       
           
       x0 = x1
       
    
   
    return[x1]

p= iterate(x0,100,Jinv1)
print(p)

def Newton(x0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    for its in range(Nmax):
       
       Jinv = Jinv1
       F = evalF(x0)
       
       x1 = x0 - Jinv.dot(F)
       
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier, its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]

p1= Newton(x0,tol1,Nmax1)
print(p1)

def Newton(x0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    for its in range(Nmax):
       J = evalJ(x0)
       Jinv = inv(J)
       F = evalF(x0)
       
       x1 = x0 - Jinv.dot(F)
       
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier, its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]

def f(x, y, z):
    return x**2 + 4*y**2 + 4*z**2 - 16

# Gradient of the function
def grad_f(x, y, z):
    return np.array([2*x, 8*y, 8*z])

# Newton's iteration scheme
def newton_iteration(x0, y0, z0, tol=1e-10, max_iter=100):
    x, y, z = x0, y0, z0
    for i in range(max_iter):
        f_val = f(x, y, z)
        grad_val = grad_f(x, y, z)
        norm_grad = np.linalg.norm(grad_val)**2
        
        if norm_grad < tol:  # Avoid division by zero in case of very small gradient
            break
        
        # Update step
        d = f_val / norm_grad
        x = x - d * grad_val[0]
        y = y - d * grad_val[1]
        z = z - d * grad_val[2]
        
        # Check for convergence
        if np.abs(f_val) < tol:
            break
    
    return x, y, z, i, f_val
print(newton_iteration(1, 1, 1, tol=1e-10, max_iter=100))


def newton_iteration_with_errors(x0, y0, z0, tol=1e-10, max_iter=100):
    x, y, z = x0, y0, z0
    errors = []
    true_point = np.array([1.093642317388195, 1.3603283832230444, 1.3603283832230444])  # Using the converged point as reference
    for i in range(max_iter):
        f_val = f(x, y, z)
        grad_val = grad_f(x, y, z)
        norm_grad = np.linalg.norm(grad_val)**2
        
        if norm_grad < tol:  # Avoid division by zero in case of very small gradient
            break
        
        # Store the current error (distance to the final converged solution)
        current_error = np.linalg.norm([x, y, z] - true_point)
        errors.append(current_error)
        
        # Update step
        d = f_val / norm_grad
        x = x - d * grad_val[0]
        y = y - d * grad_val[1]
        z = z - d * grad_val[2]
        
        # Check for convergence
        if np.abs(f_val) < tol:
            break
    
    # Return the errors across iterations
    return errors

# Perform the iteration and capture errors
errors = newton_iteration_with_errors(1, 1, 1)

# Calculate the ratio between successive errors to check quadratic convergence
error_ratios = [(errors[i+1] / (errors[i] ** 2)) for i in range(len(errors)-1)]
errors, error_ratios

print(errors, error_ratios)