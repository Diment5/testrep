import numpy as np                              #numpy
import time                                     #time (to measure the time it takes to execute code)
import math                                     #math functions
from scipy import io, integrate, linalg, signal #scipy has some linear algebra and specialized libraries
import matplotlib.pyplot as plt                 #pyplot
from matplotlib.animation import FuncAnimation  #for animation purposes
       

plt.rcParams['figure.figsize'] = [15, 5]; #changes the size of the display window

#bisection 
#1.
# biscetion code
def bisect_method(f,a,b,tol,nmax,vrb):
    #First attempt at bisection method applied to f between a and b
    
    # Initial values for interval [an,bn], midpoint xn 
    an = a; bn=b; n=0;
    xn = (an+bn)/2;
    # Current guess is stored at rn[n]
    rn=np.array([xn]); 
    r=xn;
    ier=0; 
    
    print("\n Bisection method with nmax=%d and tol=%1.1e\n" % (nmax, tol));
    
    # The code cannot work if f(a) and f(b) have the same sign. 
    # In this case, the code displays an error message, outputs empty answers and exits. 
    if f(a)*f(b)>=0:
        print("\n Interval is inadequate, f(a)*f(b)>=0. Try again \n")
        print("f(a)*f(b) = %1.1f \n" % f(a)*f(b)); 
        r = None; 
        return r
    else:
        if vrb:
            print("\n|--n--|--an--|--bn--|----xn----|-|bn-an|--|---|f(xn)|---|");
    
            fig, (ax1, ax2) = plt.subplots(1, 2)
            fig.suptitle('Bisection method results')
            ax1.set(xlabel='x',ylabel='y=f(x)')
            xl=np.linspace(a,b,100,endpoint=True); 
            yl=f(xl);
            ax1.plot(xl,yl);
    
        while n<=nmax:
            #print table row if vrb 
            if vrb:
                print("|--%d--|%1.4f|%1.4f|%1.8f|%1.8f|%1.8f|" % (n,an,bn,xn,bn-an,np.abs(f(xn))));  
            
                #################################################################################
                # Plot results of bisection on subplot 1 of 2 (horizontal). If vrb is true, pause. 
                xint = np.array([an,bn]); 
                yint=f(xint);
                ax1.plot(xint,yint,'ko',xn,f(xn),'rs');            
                #################################################################################    
                
            # If the error estimate is less than tol, get out of while loop
            if (bn-an)<2*tol: #better than np.abs(f(xn))<tol:
                #(break is an instruction that gets out of the while loop)
                ier=1; 
                
                break;  
                
            # If f(an)*f(xn)<0, pick left interval, update bn
            if f(an)*f(xn)<0:
                bn=xn;     
            else:
                #else, pick right interval, update an
                an=xn;  
       
            # update midpoint xn, increase n. 
            n += 1; 
            xn = (an+bn)/2; 
            rn = np.append(rn,xn);

    # Set root estimate to xn. 
    r=xn; 
    
    if vrb:
        ############################################################################
        # subplot 2: approximate error log-log plot
        e = np.abs(r-rn[0:n]); 
        #length of interval
        ln = (b-a)*np.exp2(-np.arange(0,e.size));
        #log-log plot error vs interval length
        ax2.plot(-np.log2(ln),np.log2(e),'r--');
        ax2.set(xlabel='-log2(bn-an)',ylabel='log2(error)');
        ############################################################################
    
    return r, rn;

#a) 
def f(x):
    x1=(x**2)*(x-1)

    return x1

tol= 10**-6
nmax=100

ba=bisect_method(f,0.5,2,tol,nmax,True)
print(ba)
#under tol=10e-6,this interval a took 20 iterations
#b)
#bb=bisect_method(f,-1,0.5,tol,nmax,True)
#this interval does not work because f(a)and f(b) are both negative
#c)
bc=bisect_method(f,-1,2,tol,nmax,True)
#this interval took 21 itertrations 

#4.2
#fixed point code
def fixed_point_method(g,x0,a,b,tol,nmax,vrb=False):
    # Fixed point iteration method applied to find the fixed point of g from starting point x0
    
    
    # Initial values for guess xn = f(xn) 
    n=0;
    xn = x0;
    # Current guess is stored at rn[n]
    rn=np.array([xn]); 
    r=xn;
    
    if vrb:
        print("\n Fixed point method with nmax=%d and tol=%1.1e\n" % (nmax, tol));
        print("\n|--n--|----xn----|---|g(xn)|---|");
        
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.suptitle('Fixed method results')
        ax1.set(xlabel='x',ylabel='y=g(x)')
        xl=np.linspace(a,b,100,endpoint=True); 
        yl=g(xl);
        ax1.plot(xl,yl); #plot of y = g(x)
        ax1.plot(xl,xl); #plot of line y = x. The intersection with g is the fixed point.
        #ax1.plot(np.array([xn,xn]),np.array([0,g(xn)]),'-rs');
    
    while n<=nmax:
        #print and pause. Remove or comment out pause if you want the code to go faster. 
        if vrb:
            print("|--%d--|%1.8f|%1.8f|" % (n,xn,np.abs(g(xn))));  
            
            #################################################################################
            # Plot results of fixed pt iteration on subplot 1 of 2 (horizontal). If vrb is true, pause. 
            ax1.plot(np.array([xn,g(xn)]),np.array([g(xn),g(xn)]),'-rs'); #horizontal line to y=x            
            #################################################################################    
                
        # If the estimate is approximately a root, get out of while loop
        if np.abs(g(xn)-xn)<tol:
            #(break is an instruction that gets out of the while loop)
            break;     
       
        # update iterate xn, increase n. 
        n += 1; 
        xn = g(xn); #apply g (fixed point step)
        
        if vrb:
            ax1.plot(np.array([xn,xn]),np.array([xn,g(xn)]),'-rs'); #vertical line back to y=g(x)
        
        rn = np.append(rn,xn); #add new guess to list of iterates

    # Set root estimate to xn. 
    r=xn; 
    
    if vrb:
        ############################################################################
        # subplot 2: approximate error log-log plot
        e = np.abs(g(rn) - rn); 
        #np.abs(r-rn[0:n]); 
        # steps array 
        ln = np.arange(0,n+1);
        #log-log plot error vs interval length
        ax2.plot(ln,np.log10(e),'r--');
        ax2.set(xlabel='n',ylabel='log10(error)');
        ############################################################################
    
    return r, rn;

def fa(x):
    return x*(1+(7-x**5)/x**2)**3
def fb(x):
    return x-(x**5-7)/(x**2-1)
def fc(x):
    return x-(x**5-7)/(5*x**4)
def fd(x):
    return x-(x**5-7)/12

(r1,r1n)=fixed_point_method(fa,1,0,8,1e-15,1000,False);
(r2,r2n)=fixed_point_method(fb,1,0,8,1e-15,1000,False);
(r3,r3n)=fixed_point_method(fc,1,0,8,1e-15,1000,False);
(r4,r4n)=fixed_point_method(fd,1,0,8,1e-15,1000,False);
e1 = np.abs(r1n-fa(r1n)); 
e2 = np.abs(r2n-fb(r2n));
e3 = np.abs(r3n-fc(r3n));
e4 = np.abs(r4n-fd(r4n))
ln1 = np.arange(0,len(e1));
ln2 = np.arange(0,len(e2));
ln3 = np.arange(0,len(e3));
ln4 = np.arrange90,len(e4);
plt.plot(ln1,np.log10(e1),'r--',label="a)");
plt.plot(ln2,np.log10(e2),'g--',label="b)");
plt.plot(ln3,np.log10(e3),'b--',label="c)");
plt.plot(ln3,np.log10(e4),'y--',label="d)");
plt.xlabel('n'); plt.ylabel('log10(error)');
plt.legend(); 

#a) and b) both have the problem of divideing by 0 if p0 is 1,
#and c),d)doenst iterate with the code