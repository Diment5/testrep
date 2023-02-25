import numpy as np
from scipy import special
import matplotlib.pyplot as plt

ti=20.
ts=-15.
a=0.138e-6
t=24*(60.**3)
def f1(x):
    return special.erf(x/(2*np.sqrt(a*t)))*(ti-ts)+ts
x1=np.linspace(0,10,100)
plt.plot(x1,f1(x1))
plt.show()

#4.b
#bisection from class
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

tol=10e-13
nmax=100

print(bisect_method(f1,0,2,tol,nmax,vrb=True))

#4.c
def newton_method(f,df,x0,tol,nmax,verb=False):
    #newton method to find root of f in [a,b]
    
    #Initialize iterates and iterate list
    xn=x0; 
    rn=np.array([x0]);
    fn=f(xn); dfn=df(xn);
    nfun=2; 
    dtol=1e-10; 
    
    if abs(dfn)<dtol:
        #If derivative is too small, Newton will fail. Error message is 
        #displayed and code terminates.
        if verb:
            print('\n derivative at initial guess is near 0, try different pt \n'); 
    else:
        n=0; 
        if verb:
            print("\n|--n--|----xn----|---|f(xn)|---|");
            
        #Iteration runs until f(xn) is small enough or nmax iterations are computed.
        
        while abs(fn)>tol and n<=nmax:
            if verb:
                print("|--%d--|%1.8f|%1.8f|" %(n,xn,np.abs(fn)));
                      
            #Compute Newton step 
            xn = xn - fn/dfn; 
            
            n+=1; 
            rn=np.append(rn,xn);
            dfn=df(xn); 
            fn=f(xn); 
            nfun+=2; 
        
        r=xn; 
        
        if abs(fn)>tol:
            print("Newton method failed to converge, n=%d, f(r)=%1.1e\n'" %(n,np.abs(fn)));
        else:
            print("Newton method converged succesfully, n=%d, f(r)=%1.1e" %(n,np.abs(fn)));
        
    return (r,rn,nfun)

def df1(x):
    return (ti-ts)*(1/(2*np.sqrt(a*t))*(2/np.sqrt(np.pi))*np.exp(-x**2/4*a*t))
#if x0=0.01
print(newton_method(f1,df1,0.01,tol,nmax,verb=False))

#if x0=2
print(newton_method(f1,df1,2,tol,nmax,verb=False))
