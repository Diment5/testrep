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