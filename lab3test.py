from mypkg.Iteration1D import Iteration1D
f = lambda x: (x-1)*(x**2)
find = Iteration1D(f,'bisection')
find.a = 0.5; find.b = 2.0
find.tol = 1e-6; find.Nmax = 100
x_bisection = find.root()
print(x_bisection)

find2 = Iteration1D(f,'bisection')
find2.a = -1.0; find2.b = 0.5
find2.tol = 1e-6; find2.Nmax = 100
x2_bisection = find2.root()
print(x2_bisection)

find3 = Iteration1D(f,'bisection')
find3.a = -1.0; find3.b = 2.0
find3.tol = 1e-6; find3.Nmax = 100
x3_bisection = find3.root()
print(x3_bisection)

find.method = 'fixedpt'
# recast problem
find.f = lambda x: (x-1)*(x**2)
# set initial guess
find.p0 = 0.8
# find the root
x_fixedpt = find.root()
print(x_fixedpt)
