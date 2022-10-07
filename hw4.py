from scipy.special import erf
import numpy as np
import matplotlib.pyplot as plt
from mypkg.Iteration1D import Iteration1D
alpha=0.138*10**-6
Ti =  20
Ts = -15
t=60**3 * 24
x = np.linspace(0,1,100)
y = erf(x/(2*(alpha*t)**1/2))*(Ti-Ts)+Ts 
plt.plot(x,y)
#plt.show()
f= lambda x:erf(x/(2*(alpha*t)**1/2))*(Ti-Ts)+Ts  
find3 = Iteration1D(f,'bisection')
find3.a = 0; find3.b = 1
find3.tol = 1e-13; find3.Nmax = 100
x3_bisection = find3.root()
print(x3_bisection)

find2 = Iteration1D(f,'newton')
find2.p0 = 1
find2.fp = lambda x: 23.35*np.exp(-0.35*x**2)
find2.tol = 1e-13; find2.Nmax = 100
x2_newton = find2.root()
print(x2_newton)