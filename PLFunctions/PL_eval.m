function[pl] = PL_eval(x,ns,s,alpha,scale)

pl_parameters;
addpath ../opoly

x = x(:);
r = x_to_r(x,scale);
beta = 2*s-2;

pl = jacobipolyn(r,ns,alpha,beta)/sqrt(scale);

rmpath ../opoly
