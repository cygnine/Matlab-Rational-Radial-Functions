function[pl] = Phi_eval(x,ns,s,alpha,scale)

phi_parameters;
addpath ../opoly

x = x(:);
r = x_to_r(x,scale);
beta = s-3/2;

pl = jacobipolyn(r,ns,alpha,beta)/sqrt(scale);

rmpath ../opoly
