import matplotlib.pyplot as plt
import numpy as np
import math
class my2DPlot:
    def __init__(self,f,a,b):
        x = np.arange(a,b,0.01)
        y = f(x)
        self.p = plt.plot(x,y)
    def show(self):
        plt.show()
    def dotted(self):
        self.p[0].set_linestyle('dotted')
    def labels(self,x,y):
        plt.xlabel(x)
        plt.ylabel(y)
    def color(self,colorName):
        self.p[-1].set_color(colorName)
    def logy(self):
        plt.yscale('log')
    def logx(self):
        plt.xscale('log')
    def savefig(self,fileName):
        plt.savefig(fileName)

