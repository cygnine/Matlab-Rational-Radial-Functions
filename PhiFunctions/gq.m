function[x,w] = gq(N,s,alpha,scale)

phi_parameters;
addpath ../opoly
beta = s-3/2;

[r,w] = jacobipoly_gq(N,alpha,beta);

x = r_to_x(flipud(r),scale);
w = flipud(w)*scale;

rmpath ../opoly
