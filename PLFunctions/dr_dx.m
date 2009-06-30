function[drdx] = dr_dx(x,scale)

pl_parameters;
standard_scale;

drdx = 1/scale*-2./(1+x).^2;
