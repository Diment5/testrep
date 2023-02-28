import numpy as np
from numpy import random as rand
import matplotlib.pyplot as plt

#1 G'(X)=f(x)*f''(x)/(f'(x))**2 < 1

#2
def combine(f,df,ddf,a,b,tol,Nmax):
    '''
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
    '''

    '''     first verify there is a root we can find in the interval '''
    p = np.zeros(Nmax+1);

    fa = f(a); fb = f(b);
    g=lambda x: f(x)*ddf(x)/(df(x))**2
    ######################################
    # if conditions for bisection aren't met, exit (ier=1)
    if (fa*fb>0):
        ier = 1;
        astar = a;
        p[0] = astar;
        return [p,astar, ier,count]

    ######################################
    # if f(a)=0 or f(b)=0, exit (ier=0)
    ''' verify end point is not a root '''
    if (fa == 0):
        astar = a;
        ier =0;
        p[0] = astar;
        return [p,astar, ier,count]

    if (fb ==0):
        astar = b;
        ier = 0;
        p[0] = astar;
        return [p,astar, ier,count]
    #####################################

    # iteration starts (while count<Nmax)
    count = 0;
    while (count < Nmax):
        c = 0.5*(a+b);
        fc = f(c);
        p[count+1] = c;

        # if midpoint is a zero, exit (ier=0)
        if (fc ==0):
            astar = c;
            ier = 0;
            p[count+1] = astar;
            return [p,astar, ier,count]

        if (fa*fc<0):
            b = c;
        elif (fb*fc<0):
            a = c;
            fa = fc;
        else:
            # if fc=0, exit (ier=3)
            astar = c;
            ier = 3;
            p[count+1] = astar;
            return [p,astar, ier,count]

        # once our interval is smaller than tol, exit (ier=0)
        if (abs(g(c))< 1.):# our basin step and appending the newton method
            pn = np.zeros(Nmax+1);
            pn[0] = c
            for it in range(Nmax):
                p1 = c-f(c)/df(c)
                p[it+1] = p1
                if (abs(p1-c) < tol):
                    pstar = p1
                    info = 0
                    return [p,pstar,info,it,"the combine method gives us the root of",pstar]
                c = p1;
                pstar = p1;
                info = 1;
            return [p,pstar,info,it,"bisection fail"]

            
        
    
        
        if (abs(b-a)<tol):
            astar = a;
            ier =0;
            p[count+1] = astar;
            return [p,astar, ier,count]

        count = count +1;

    # If we are here, our algorithm ran out of iter. exit (ier=2)
    astar = a;
    ier = 2;
    return [p,astar,ier,count];
    

    


#3, we also need to input df and ddf

#5 this method is garentee to converge and is faster than the regular bisection.
#its limit is that it needs the user to input the frist and second order derivative

#6 
def f6(x):
    return np.exp(x**2+7*x-30)-1
tol=10e-13
Nmax=100

#6.A bisection 
def bisection(f,a,b,tol,Nmax):
    '''
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
    '''

    '''     first verify there is a root we can find in the interval '''
    p = np.zeros(Nmax+1);

    fa = f(a); fb = f(b);
    ######################################
    # if conditions for bisection aren't met, exit (ier=1)
    if (fa*fb>0):
        ier = 1;
        astar = a;
        p[0] = astar;
        return [p,astar, ier,count]

    ######################################
    # if f(a)=0 or f(b)=0, exit (ier=0)
    ''' verify end point is not a root '''
    if (fa == 0):
        astar = a;
        ier =0;
        p[0] = astar;
        return [p,astar, ier,count]

    if (fb ==0):
        astar = b;
        ier = 0;
        p[0] = astar;
        return [p,astar, ier,count]
    #####################################

    # iteration starts (while count<Nmax)
    count = 0;
    while (count < Nmax):
        c = 0.5*(a+b);
        fc = f(c);
        p[count+1] = c;

        # if midpoint is a zero, exit (ier=0)
        if (fc ==0):
            astar = c;
            ier = 0;
            p[count+1] = astar;
            return [p,astar, ier,count]

        if (fa*fc<0):
            b = c;
        elif (fb*fc<0):
            a = c;
            fa = fc;
        else:
            # if fc=0, exit (ier=3)
            astar = c;
            ier = 3;
            p[count+1] = astar;
            return [p,astar, ier,count]

        # once our interval is smaller than tol, exit (ier=0)
        if (abs(b-a)<tol):
            astar = a;
            ier =0;
            p[count+1] = astar;
            return [p,astar, ier,count]

        count = count +1;

    # If we are here, our algorithm ran out of iter. exit (ier=2)
    astar = a;
    ier = 2;
    return [p,astar,ier,count];
 
print(bisection(f6,2,4.5,tol,Nmax))
#the bisection took 41 iterations 
#the cost of each itertion is to evaluate f(x) twice and divided it by two

#6.B newton
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
    p = np.zeros(Nmax+1);
    p[0] = p0
    for it in range(Nmax):
        p1 = p0-f(p0)/fp(p0)
        p[it+1] = p1
        if (abs(p1-p0) < tol):
            pstar = p1
            info = 0
            return [p,pstar,info,it]
        p0 = p1;
        pstar = p1;
        info = 1;
    return [p,pstar,info,it]

df6= lambda x: (2*x+7)*np.exp(x**2+7*x-30)
print(newton(f6,df6,4.5,tol, Nmax))
# the newton method took 26 iterations
# the cost of each iterations is evaluating f and df 

#6.C combine 
ddf6 = lambda x: 2*(2*x+7)*np.exp(x**2+7*x-30)+(2*x+7)**2*np.exp(x**2+7*x-30)
print(combine(f6,df6,ddf6,2,4.5,tol,Nmax))
#the combine method took 6 ietrations, for each itertions the cost is the sum of both bisection 
#and newton. It converges the faestest, and is the most cost effective.