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