classdef my2DPlot < handle
properties
p
end
methods
function self = my2DPlot(f,a,b)
x = linspace(a,b,100);
y = f(x);
self.p = [plot(x,y)];
end
function [] = dotted(self)
self.p(1).LineStyle = 'none';
self.p(1).Marker = '.';
end
end
