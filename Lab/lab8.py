import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.special import comb
from scipy.interpolate import splev, splrep, splint,CubicSpline; 
plt.rcParams['figure.figsize'] = [15, 8];

def driver():


    f = lambda x: np.exp(x)

    N = 10;
    ''' interval'''
    a = 0;
    b = 1;

    ''' create equispaced interpolation nodes'''
    xint = np.linspace(a,b,N+1);

    ''' create interpolation data'''
    yint = f(xint);

    Neval = 1000;
    xeval = np.linspace(a,b,Neval+1);

    ''' Linear spline evaluation '''
    yeval_ls = eval_lin_spline(xeval,xint,yint,N);

    ''' create points for evaluating the Lagrange interpolating polynomial'''
    yeval_l= np.zeros(Neval+1)
    yeval_dd = np.zeros(Neval+1)

    '''Initialize and populate the first columns of the
     divided difference matrix. We will pass the x vector'''
    y = np.zeros( (N+1, N+1) )

    for j in range(N+1):
       y[j][0]  = yint[j]

    y = dividedDiffTable(xint, y, N+1)
    ''' evaluate lagrange poly '''
    for kk in range(Neval+1):
       yeval_l[kk] = eval_lagrange(xeval[kk],xint,yint,N)
       yeval_dd[kk] = evalDDpoly(xeval[kk],xint,y,N)

    ''' create vector with exact values'''
    fex = f(xeval)

    plt.figure()
    plt.plot(xeval,fex,'ro-')
    plt.plot(xeval,yeval_l,'bs--')
    plt.plot(xeval,yeval_dd,'c.--')
    plt.plot(xeval,yeval_ls,'g--')
    plt.legend()

    plt.figure()
    err_l = abs(yeval_l-fex)
    err_dd = abs(yeval_dd-fex)
    err_ls = abs(yeval_ls-fex)
    plt.semilogy(xeval,err_l,'ro--',label='lagrange')
    plt.semilogy(xeval,err_dd,'bs--',label='Newton DD')
    plt.semilogy(xeval,err_ls,'g--',label='lin spline')
    plt.legend()
    plt.show()

def eval_line(x,x0,y0,x1,y1):
    lin = (1/(x1-x0))*(y0*(x1-x) + y1*(x-x0));
    return lin;

def find_int(xeval,a,b):
    ind = np.where(np.logical_and(xeval>=a,xeval<=b));
    return ind;

def eval_lin_spline(xeval,xint,yint,N):
    Neval = len(xeval);
    yeval = np.zeros(Neval);

    for n in range(N):
        indn = find_int(xeval,xint[n],xint[n+1]);
        yeval[indn] = eval_line(xeval[indn],xint[n],yint[n],xint[n+1],yint[n+1]);

    return yeval;

def eval_lagrange(xeval,xint,yint,N):

    lj = np.ones(N+1)

    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.

    for jj in range(N+1):
       yeval = yeval + yint[jj]*lj[jj]

    return(yeval)


''' create divided difference matrix'''
def dividedDiffTable(x, y, n):

    for i in range(1, n):
        for j in range(n - i):
            y[j][i] = ((y[j][i - 1] - y[j + 1][i - 1]) /
                                     (x[j] - x[i + j]));
    return y;

def evalDDpoly(xval, xint,y,N):
    ''' evaluate the polynomial terms'''
    ptmp = np.zeros(N+1)

    ptmp[0] = 1.
    for j in range(N):
      ptmp[j+1] = ptmp[j]*(xval-xint[j])

    '''evaluate the divided difference polynomial'''
    yeval = 0.
    for j in range(N+1):
       yeval = yeval + y[0][j]*ptmp[j]

    return yeval

#3.2
def eval_lin_spline(xeval,xint,yint,N):
    Neval = len(xeval);
    yeval = np.zeros(Neval);

    for n in range(N):
        indn = find_int(xeval,xint[n],xint[n+1]);
        yeval[indn] = eval_line(xeval[indn],xint[n],yint[n],xint[n+1],yint[n+1]);

    return yeval;
def find_int(xeval,a,b):
    ind = np.where(np.logical_and(xeval>=a,xeval<=b));
    return ind;

def eval_line(x,x0,y0,x1,y1):
    lin = (1/(x1-x0))*(y0*(x1-x) + y1*(x-x0));
    return lin;

f1=lambda x: 1/(1+(10*x)**2)
Neval = 1000
N= 10
xeval = np.linspace(-1,1,Neval+1)
xint = np.linspace(-1,1,N+1)
yint = f1(xint)
print(eval_lin_spline(xeval,xint,yint,N))

fex = f1(xeval)
yeval_ls = eval_lin_spline(xeval,xint,yint,N);

plt.figure()
plt.plot(xeval,fex,'ro-')
#plt.plot(xeval,yeval_l,'bs--')
#plt.plot(xeval,yeval_dd,'c.--')
plt.plot(xeval,yeval_ls,'g--')
plt.legend()

plt.figure()
#err_l = abs(yeval_l-fex)
#err_dd = abs(yeval_dd-fex)
err_ls = abs(yeval_ls-fex)
#plt.semilogy(xeval,err_l,'ro--',label='lagrange')
#plt.semilogy(xeval,err_dd,'bs--',label='Newton DD')
plt.semilogy(xeval,err_ls,'g--',label='lin spline')
plt.legend()
plt.show()

#This behaves more accurate and the error is even across each interval 

#3.2

def g(x):
    return 1/(1+10*x**2);
def dg(x):
    return -20*x/(1+10*x**2)**2;
def dg2(x):
    return -20/(1+10*x**2)**2 + (800*x**2)/(1+10*x**2)**3
def dg3(x):
    return -20*dg(x) + (1600*x)/(1+10*x**2)**3 + (800*x**2)*(-3*20*x)/(1+10*x**2)**4;
def cubicspline(fun,n,der=0,dfun=None):
    # Interpolation pts are n+1 equispaced
    xi=np.linspace(-1,1,n+1);
    h = 2/n;     # uniform spacing h
    f = fun(xi); # function values
    # Initialize two subplots  
    fig, (ax1, ax2)  = plt.subplots(1, 2);
    
    
    # Initialize matrix A and right-hand-side rhs
    A = np.zeros((n-1,n-1)); #A is (n-1) x (n-1)
    rhs = np.zeros(n-1);     #rhs is (n-1) x 1
    ct = (6/(h**2)); #This constant appears in the rhs
    
    for i in np.arange(n-1):
        # (6/h^2) * (f[i] - 2f[i+1] + f[i+2])
        rhs[i] = ct*(f[i]-2*f[i+1]+f[i+2]);
        if i == 0:
        # Matrix A has 3 diagonals with values [1,4,1]
            # special case for first row [4 1 0 ... 0]
            A[0,0:2]=[4,1];
        elif i==n-2:
            # special case for last row [0 0 ... 1 4]
            A[n-2,n-3:n-1]=[1,4];
        else:
            # general case [0 0 ... 1 4 1 ... 0 0]
            A[i,i-1:i+2]=[1,4,1];
    
    # We create coefficient vector a of size n+1
    a = np.zeros(n+1); 
    # (a_1,...,a_n-1) are given by the solution A^{-1}*rhs. 
    # a0 = an = 0 by the natural spline condition. 
    tmp = np.linalg.solve(A,rhs); 
    a[1:n] = tmp; 
    # b and c coefficients are functions of f(xi) and coeffs a. 
    b = f[0:n] - (1/ct)*a[0:n]; #bi = fi - (h^2/6)*ai
    c = f[1:n+1] - (1/ct)*a[1:n+1]; #ci = f_{i+1} - (h^2/6)*a_{i+1}
    for i in np.arange(1,n+1):
        y = np.linspace(xi[i-1],xi[i],100);
        
            # Original spline formula
        si = (a[i-1]/(6*h))*(xi[i]-y)**3 + (a[i]/(6*h))*(y-xi[i-1])**3 \
            + (b[i-1]/h)*(xi[i]-y) + (c[i-1]/h)*(y-xi[i-1]);
        return si
print(cubicspline(f1,10,der=0,dfun=None))
def eval_cub_spline(xeval,xint,yint,N):
    Neval = len(xeval);
    yeval = np.zeros(Neval);

    for n in range(N):
        indn = find_int(xeval,xint[n],xint[n+1]);
        yeval[indn] = cubicspline(f1,10,der=0,dfun=None)[n]

    return yeval;


print(eval_cub_spline(xeval,xint,yint,N))
xeval1 = cubicspline(f1,10,der=0,dfun=None)


    