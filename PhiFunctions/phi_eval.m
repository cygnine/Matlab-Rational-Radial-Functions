function[pl] = phi_eval(x,ns,s,alpha,scale)

phi_parameters;
addpath ../opoly
warning('off','MATLAB:rmpath:DirNotFound');

x = x(:);
N = size(x,1);
r = x_to_r(x,scale);
beta = s-3/2;

%pl = jacobipolyn(r,ns,alpha,beta)/sqrt(scale);
pl = Phi_eval(x,ns,s,alpha,scale);
weight = (1+r).^(s/2);

pl = spdiags(weight,0,N,N)*pl;

rmpath ../opoly
warning('on','MATLAB:rmpath:DirNotFound');
