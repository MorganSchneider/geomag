function [x0] = interpolate(dF1,dF2,x1,x2)

x0 = x1 - dF1 * ((x2 - x1) / (dF2 - dF1));

return