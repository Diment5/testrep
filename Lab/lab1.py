import numpy as np
import matplotlib.pyplot as plt

#exercise 3.2
#1
x = np.linspace(0,10,11)
y = np.array([0,1,2,3,4,5,6,7,8,9,10])
#2,3
print("the first three entires of x are", x[0:3])
#4,5
w = 10**(-np.linspace(1,10,10))
s=3*w
x1 = np.linspace(1,10,10)
plt.semilogy(x1,w)
plt.semilogy(x1,s)
plt.xlabel("x1")
plt.ylabel("w,s")
plt.show()

#excercise 4.2
#1, with samplecodes
import numpy as np
import numpy.linalg as la
import math
def driver():
    n = 100
    x = np.linspace(0,np.pi,n)
# this is a function handle. You can use it to define
# functions instead of using a subroutine like you
# have to in a true low level language.
    f = lambda x: x**2 + 4*x + 2*np.exp(x)
    g = lambda x: 6*x**3 + 2*np.sin(x)
    y = f(x)
    w = g(x)
# evaluate the dot product of y and w
    dp = dotProduct(y,w,n)
# print the output
    print('the dot product is : ', dp)
    
    return

def dotProduct(x,y,n):
    dp = 0.
    for j in range(n):
        dp = dp + x[j]*y[j]
    
    return dp

driver()
 # i am not sure of what the question is aking about change the vectors
 # in the dot product code to ones that are orthogonal.

#2
def mmulti(x,y): 
    xrow=len(x)  #find matirx parameters
    xcol=len(x[0])
    yrow=len(y)
    ycol=len(y[0])
    mp=np.zeros((xrow,ycol)) #set result matrix
    for i in range(xrow): #loop through all the colum and rows
        for j in range(xcol):
            for k in range(ycol):
                mp[i][k] += x[i][j]*y[j][k]
    
    return mp
#test with 2 2by2 matrix (correct)
x1=np.array([[1,2],[3,4]])
y1=np.array([[5,6],[7,8]])
m1=mmulti(x1,y1)
print(m1)
 
#test with bigger 3by4 and 4by5 matirx (correct)
x2=np.array([[8,7,4,2],[9,2,5,1],[2,5,6,2]])
y2=np.array([[9,5,6,7,3],[2,3,6,4,5],[3,6,7,4,7],[3,5,6,7,8]])
m2=mmulti(x2,y2)
print(m2)

#3 matrix multipy in numpy
nm1=np.matmul(x1,y1)
nm2=np.matmul(x2, y2)

print(nm1)
print(nm2)
# although i dont know which runs fatser
# presumabaly still numpy
# i feel like an idiot looping through all the index 
# after seeing that one function from numpy
#testrep2 push

