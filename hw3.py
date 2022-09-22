from itertools import count
from mypkg.Iteration1D import Iteration1D
import numpy as np
from mypkg.my2DPlot import my2DPlot
f = lambda x: x**3 +x -4
find = Iteration1D(f,'bisection')
find.a = 1; find.b = 3
find.tol = 1e-3; find.Nmax = 100
x_bisection = find.root()
x_iteration = find.count
print(x_bisection,x_iteration)




find.method = 'fixedpt'
# need to specify initial guess for this method
find.p0 = 2.3
# recasted problem
find.f = lambda x: -np.sin(2*x)-3/4 +5*x/4
find.tol = 1e-10; find.Nmax = 100
# find root
x_fixedpt = find.root()
print(x_fixedpt)
