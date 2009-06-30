function[x,w] = grq(N,s,alpha,scale)

pl_parameters;
addpath ../opoly
beta = 2*s-2;

[r,w] = jacobipoly_grq(N,alpha,beta,1);

x = r_to_x(flipud(r),scale);
w = flipud(w)*scale;

rmpath ../opoly
