import numpy as np
import matplotlib.pyplot as plt


def compute_order(x,xstar):

    diff1=np.abs(x[1::]-xstar)
    diff2=np.abs(x[0:-1]-xstar)
    fit=np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()),1)
    print('the order equation is')
    print('log(|p_{n+1}-p|)=log(lamdba)+alpha*log(|p_n-p|)')
    print('lambda=', str(np.exp(fit[1])))
    print('alpha=', str(fit[0]))

    return[fit, diff1, diff2]



#prelab
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    # make an array of zeros of length Nmax
    x = np.zeros((Nmax,1))
    # save the initial guess
    x[0] = x0
    while (count < Nmax):
       count = count +1
       x1 = f(x0)
       # save the current iterate
       x[count] = x1
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          # truncate the array to have only count entries
          x = x[0:count]
          return [xstar,x,ier,count]
       x0 = x1

    xstar = x1
    x = x[0:count]
    ier = 1
    return [xstar,x,ier,count]
def g(x):
    return np.sqrt(10/(x+4))
(rstar1,r1,er1,rn1)=fixedpt(g,1.5,1e-10,100)
print(r1,rn1)
# it take 15 iterations 


#order of convergence 
od1= compute_order(r1[0:-1],r1[-1])
print(od1)


#3.1
#p = ((-pn1**2+pn*pn2)/(pn-2*pn1+pn2))
#3.1
def px(x):
    phat=np.zeros(len(x))
    i=0
    while i < len(x)-2:
        phat[i]+=x[i]-((x[i+1]-x[i])**2/(x[i]-2*x[i+1]+x[i+2]))
        i += 1
        
    return phat

print(px(r1))
p0= 1.3652300134140976

#the fixpt took 12 iterations, wherease the aitken's took only 3

#3.4
#define stf(a,b,c):a-(b-a)**2/(c-2b-a)
#input: p0, g(x),nmax, and tol

#1. count=1, 
#2. x= array of zeros of size nmax,x[0]=p0,
#3. while count < nmax
#4. a=x[count-1],b=g(a),c=g(b)
#5. x[count] = stf(a,b,c)
#6. return x
#3.3.2

def stf(a,b,c):
    return a-(b-a)**2/(c-2*b-a)

def stfn(g,p0,tol,nmax):
    x=np.zeros(nmax)
    count=1
    x[0]=p0
    while count < nmax and np.abs(g(x[count])-x[count])<tol:
        a=x[count-1]
        b=g(a)
        c=g(b)
        x[count]=stf(a,b,c)
        count+= 1 
    
    
    return x
    
px2=stfn(g,p0,1e-10,100)
print(px2)
plt.semilogy(range(len(r1)),r1-p0)
plt.semilogy(range(len(px(r1))),px(r1)-p0)
plt.semilogy(range(len(px2),px2-p0))
plt.show()
# the latter two mathob converge at the same rate, they are both faster than the first fixedpt