function[pl] = dpl_eval(x,ns,s,alpha,scale)

pl_parameters;
x = x(:);
N = size(x,1);
r = x_to_r(x,scale);

pl = PL_eval(x,ns,s,alpha,scale);
dpl = dPL_eval(x,ns,s,alpha,scale);
drdx = dr_dx(x,scale);
weight = 1/sqrt(2)*(1+r).^s;
dweight = 1/sqrt(2)*s*(1+r).^(s-1).*drdx;

% d(w*PL) = d(PL)*w + PL*d(w)
pl = spdiags(weight,0,N,N)*dpl + spdiags(dweight,0,N,N)*pl;
