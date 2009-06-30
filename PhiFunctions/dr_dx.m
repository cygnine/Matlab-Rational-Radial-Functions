function[drdx] = dr_dx(x,scale)

phi_parameters;
standard_scale;

drdx = -4*x./(1+x.^2).^2/scale;
