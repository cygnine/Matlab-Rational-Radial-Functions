function[pl] = pl_eval(x,ns,s,alpha,scale)

pl_parameters;
addpath ../opoly
warning('off','MATLAB:rmpath:DirNotFound');

x = x(:);
N = size(x,1);
r = x_to_r(x,scale);
beta = 2*s-2;

%pl = jacobipolyn(r,ns,alpha,beta)/sqrt(scale);
pl = PL_eval(x,ns,s,alpha,scale);
weight = 1/sqrt(2)*(1+r).^s;

pl = spdiags(weight,0,N,N)*pl;

rmpath ../opoly
warning('on','MATLAB:rmpath:DirNotFound');
