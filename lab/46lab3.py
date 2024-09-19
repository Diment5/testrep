# import libraries
import numpy as np



# use routines    
 #   f = lambda x: x**3+x-4
  #  a = 1
  #  b = 4

#    f = lambda x: np.sin(x)
#    a = 0.1
#    b = np.pi+0.1

tol = 1e-6

    #[astar,ier] = bisection(f,a,b,tol)
  #  print('the approximate root is',astar)
   # print('the error message reads:',ier)
   # print('f(astar) =', f(astar))




# define routines
def bisection(f,a,b,tol):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 

    fa = f(a)
    fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier]

    count = 0
    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    return [astar, ier]
      
#now lab3 
#4.1
#a)


f1= lambda x: (x**2)*(x-1)
a1=0.5
b1=2.
[astar1,ier1] = bisection(f1,a1,b1,tol)
print('the approximate root for a is',astar1)
print('the error message reads:',ier1)
print('f1(astar1) =', f1(astar1))

#b)

a2=-1.
b2=0.5
[astar2,ier2] = bisection(f1,a2,b2,tol)
print('the approximate root for b is',astar2)
print('the error message reads:',ier2)
print('f1(astar2) =', f1(astar2))

#c)

a3=-1
b3=2.
[astar3,ier3] = bisection(f1,a3,b3,tol)
print('the approximate root for a is',astar3)
print('the error message reads:',ier3)
print('f1(astar3) =', f1(astar3))

# a) and c) works, but b) doesnt, this is becasue f(-1)=0, which returned error code 0
#it is possi

#4.2
#a)
f2a=lambda x: (x-1)*(x-3)*(x-5)
a21=0
b21=2.4
tol2=1e-5
[astar21,ier21] = bisection(f2a,a21,b21,tol2)
print('the approximate root for 2a is',astar21)
print('the error message reads:',ier21)
print('f1(astar21) =', f2a(astar21))

#b)
f2b=lambda x: ((x-1)**2)*(x-3)
a22=0
b22=2.
tol1=1e-5
[astar22,ier22] = bisection(f2b,a22,b22,tol1)
print('the approximate root for 2b is',astar22)
print('the error message reads:',ier22)
print('f1(astar22) =', f2b(astar22))

#c)
f2c=lambda x: np.sin(x)
a23=0
b23=.1
tol1=1e-5
[astar23,ier23] = bisection(f2c,a23,b23,tol1)
print('the approximate root for 2c is',astar23)
print('the error message reads:',ier23)
print('f1(astar22) =', f2c(astar23))