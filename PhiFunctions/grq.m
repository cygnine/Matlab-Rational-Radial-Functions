function[x,w] = grq(N,s,alpha,scale)

phi_parameters;
addpath ../opoly
beta = s-3/2;

[r,w] = jacobipoly_grq(N,alpha,beta,1);

x = r_to_x(flipud(r),scale);
w = flipud(w)*scale;

rmpath ../opoly
