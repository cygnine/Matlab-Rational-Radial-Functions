function[drdx] = dr_dx_divided_x(x,scale)

phi_parameters;
standard_scale;

drdx = -4./(1+x.^2).^2/scale^2;

% The extra 1/scale is because dividing by x needs to be scaled appropriately
