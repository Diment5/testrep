#prelab
import numpy as np
def  eval_legendre(n,x):
    a=np.zeros(n+1)
    a[0]=1
    a[1]=x
    for j in range(2,n):
        a[j]=1/(1+n)*((2*n+1)*x*a[j-1]-j*a[j-2])
        j += 1
    return a 

#labday
from scipy.integrate import quad
def aj(f,w,n,x,a,b):
    phi=np.zeros(n+1)
    phi[0]=1
    phi[1]=x
    for j in range(2,n):
        phi[j]=1/(1+n)*((2*n+1)*x*a[j-1]-j*a[j-2])
        j += 1
    r1=np.zeros(n+1)
    r2=np.zeros(n+1)
    a=np.zeros(n+1)
    for i in range(0,n):
        r1[i]=quad(phi[i]*f(x)*w(x),a,b)
        r2[i]=quad((phi[i]**2.)*w(x))
        
        i+=1
    a=r1/r2
    return a

w=1
a=-1
b=1
def eval_legendre_expansion(f,a,b,w,n,x):
    p=np.zeros(n)
    for j in range(0,n):
        p[j]=eval_legendre(n,x)[j]*aj(f,w,n,x,a,b)[j]
        j+=1
    return sum(p)
f = lambda x: 1/(1+x**2)




