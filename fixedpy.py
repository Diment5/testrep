from mypkg.Iteration1D import Iteration1D
f = lambda x: (x-1)*(x**2)
find = Iteration1D(f,'fixedpt')
find.method = 'fixedpt'
# need to specify initial guess for this method
find.p0 = 1.2 
# recasted problem
find.tol = 1e-6; find.Nmax = 100
# find root
x_fixedpt = find.root()
print(x_fixedpt)