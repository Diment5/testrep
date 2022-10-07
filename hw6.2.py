import mypkg.prini as prini
import numpy as np
import math
import time

from numpy.linalg import inv 
from numpy.linalg import norm 
x01 = np.array([0., 0., 0.])
x02 = np.array([1., 1., 1.])
    
Nmax1 = 100
tol1 = 1e-10
    
     
def evalF(x): 
    F = np.zeros(3)
    
    F[0] = x[0]+np.cos(x[0]*x[1]*x[2])-1
    F[1] = (1-x[0])**(1/4)+x[1]+0.05*x[2]-0.15*x[2]-1
    F[2] = -x[0]**2-0.1*x[1]**2+0.01*x[1]+x[2]-1
    return F
    
def evalJ(x): 
    
    J = np.array([[1-x[1]*x[2]*np.sin(x[0]*x[1]*x[2]), -x[2]*x[0]*np.sin(x[0]*x[1]*x[2]), -x[0]*x[1]*np.sin(x[0]*x[1]*x[2])], 
        [1/4 *(1.-x[0])**(-3/4), 1., 0.1*x[2]-0.15], 
        [-2*x[0], -0.2*x[1]+0.01, 1.]])
    return J
def Nwton(x0,tol,Nmax):
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
 
nt=Nwton(x01,tol1,Nmax1)
print("newton",nt)

def LazyNewton(x0,tol,Nmax):
    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    J = evalJ(x0)
    Jinv = inv(J)
    for its in range(Nmax):
       F = evalF(x0)
       x1 = x0 - Jinv.dot(F)
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier,its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]   
    


lz= LazyNewton(x01,tol1,Nmax1) 
print("lazynewton",lz)
    
def Broyden(x0,tol,Nmax):
    '''tol = desired accuracy
    Nmax = max number of iterations'''
    '''Sherman-Morrison 
   (A+xy^T)^{-1} = A^{-1}-1/p*(A^{-1}xy^TA^{-1})
    where p = 1+y^TA^{-1}Ax'''
    '''In Newton
    x_k+1 = xk -(G(x_k))^{-1}*F(x_k)'''
    '''In Broyden 
    x = [F(xk)-F(xk-1)-\hat{G}_k-1(xk-xk-1)
    y = x_k-x_k-1/||x_k-x_k-1||^2'''
    ''' implemented as in equation (10.16) on page 650 of text'''
    
    '''initialize with 1 newton step'''
    
    A0 = evalJ(x0)
    v = evalF(x0)
    A = np.linalg.inv(A0)
    s = -A.dot(v)
    xk = x0+s
    for  its in range(Nmax):
       '''(save v from previous step)'''
       w = v
       ''' create new v'''
       v = evalF(xk)
       '''y_k = F(xk)-F(xk-1)'''
       y = v-w;                   
       '''-A_{k-1}^{-1}y_k'''
       z = -A.dot(y)
       ''' p = s_k^tA_{k-1}^{-1}y_k'''
       p = -np.dot(s,z)                 
       u = np.dot(s,A) 
       ''' A = A_k^{-1} via Morrison formula'''
       tmp = s+z
       tmp2 = np.outer(tmp,u)
       A = A+1./p*tmp2
       ''' -A_k^{-1}F(x_k)'''
       s = -A.dot(v)
       xk = xk+s
       if (norm(s)<tol):
          alpha = xk
          ier = 0
          return[alpha,ier,its]
    alpha = xk
    ier = 1
    return[alpha,ier,its]

byd= Broyden(x01,tol1,Nmax1)

print("broyden",byd)

def evalg(x02):
    F = evalF(x02)
    G=F[0]**2 + F[1]**2 + F[2]**2
    return G


def Steepestdescent(n,x02,tol,Nmax):
    F = evalF(x02)
    J = evalJ(x02)
    

    k=1
    while k <= Nmax:
        g1= evalg(x02)
        z=2*np.transpose(J).dot(F)
        z0=norm(z)
        if z0 := 0:
            xstar=x02
            gstar=g1
            return["zerogradient", xstar, gstar]
            break
        z = z/z0 
        a1=0
        a3=1
        g3=evalg(x02-a3*x02)
        
        while g3 >= g1:
            a3=a3/2
            g3=evalg(x02-a3*x02)
            if a3 < tol/2:
                xstar=x02
                gstar=g1
                return("not likely to imporve",xstar,gstar)
            a2=a3/2
            g2=evalg(x02-a2*x02)
            h1=(g2-g1)/a2
            h2=(g3-g2)/(a3-a2)
            h3=(h2-h1)/a3
            a0=0.5*(a2-h1/h3)
            g0=evalg(x02-a0*z)
            if g0 > g3:
                a=a3
                g=g3
            elif g0 < g3:
                a=a0
                g=g0
            x02=x02-a*x02
            if np.abs(g-g1) < tol:
                xstar=x02
                return[xstar,g]
            k=k+1

        return["maximum iteration exceeded"]
        break
    

sd= Steepestdescent(3,x02,tol1,Nmax1)
print("Steepestdescent",sd)

        
    
    

        
    



