import numpy as np
from mypkg.Iteration1D import Iteration1D

f= lambda x: np.exp(x**2-7*x+30)-1
fp = lambda x: np.exp(x**2-7*x+30)*(2*x+7)
fpp = lambda x: np.exp(x**2-7*x+30)*(4*x**2+28*x+51)
find2 = Iteration1D(f,'bisection')
find2.a = 2; find2.b = 4.5
find2.tol = 1e-8; find2.Nmax = 100
x2_bisection = find2.root()
print(x2_bisection)

find3 = Iteration1D(f,'newton')
find3.p0 = 4.5 
find3.fp = fp 
find3.tol = 1e-8; find3.Nmax = 100
x3_bisection = find3.root()
print(x3_bisection)


tol=1e-8
Nmax = 100

def bisection(f,fp,fpp,a,b,tol,Nmax):
    """
    Inputs:
      f,a,b       - function and endpoints of initial interval
      tol, Nmax   - bisection stops when interval length < tol
                  - or if Nmax iterations have occured
    Returns:
      astar - approximation of root
      ier   - error message
            - ier = 1 => cannot tell if there is a root in the interval
            - ier = 0 == success
            - ier = 2 => ran out of iterations
            - ier = 3 => other error ==== You can explain
    """

    '''
     first verify there is a root we can find in the interval
    '''
    
    
    gp = f*fpp/fp**2 -0.98
    count = 0
    ga = gp(a); gb = gp(b)
    if (fa*gb>0):
       ier = 1
       astar = a
       return [astar, count]

    '''
     verify end point is not a root
    '''
    if (ga == 0):
      astar = a
      ier =0
      return [astar, count]

    if (gb ==0):
      astar = b
      ier = 0
      return [astar, count]

    
    xb = np.zeros((Nmax,1))
    xb[0] = 0.5*(a+b)
    while (count < Nmax):
      c = 0.5*(a+b)
      gpc = gp(c)

      

      if (gpc ==0):
        astar = c
        ier = 0
        return [astar, count]

      
      elif(ga*gpc<0):
         b = c
      elif (gb*gpc<0):
        a = c
        fa = gpc
      else:
        astar = c
        ier = 3
        return [astar, count]

      if (abs(b-a)<tol):
        astar = a
        ier =0
        return [astar, count]
      
      count = count +1
      
    
    astar = a
    ier = 2
    return [astar, count]
p0= bisection[0]
print(p0)

def newton(f,fp,p0,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+1)
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

aprx = newton

print(aprx)
