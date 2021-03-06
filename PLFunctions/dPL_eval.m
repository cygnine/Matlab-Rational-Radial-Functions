function[pl] = dPL_eval(x,ns,s,alpha,scale)

pl_parameters;
addpath ../opoly

x = x(:);
beta = 2*s-2;
r = x_to_r(x,scale);
drdx = dr_dx(x,scale);

pl = djacobipolyn(r,ns,alpha,beta)/sqrt(scale);
N = size(x,1);

pl = spdiags(drdx,0,N,N)*pl;

rmpath ../opoly
