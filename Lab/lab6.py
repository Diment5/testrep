import numpy as np

#prelab
#1
f= lambda x: np.cos(x)
h = 0.01 *2. **(-np.arange(0,10))
df1= lambda x: (f(x+h)-f(x))/h
df2= lambda x: (f(x+h)-f(x-h))/(2*h)
x1= np.pi/2
print(df1(x1))
print(df2(x1))
#2 both method are linear