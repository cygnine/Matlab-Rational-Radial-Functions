function[pl] = dphi_eval(x,ns,s,alpha,scale)

phi_parameters;
x = x(:);
N = size(x,1);
r = x_to_r(x,scale);

pl = Phi_eval(x,ns,s,alpha,scale);
dpl = dPhi_eval(x,ns,s,alpha,scale);
drdx = dr_dx(x,scale);
weight = (1+r).^(s/2);
dweight = (s/2)*(1+r).^(s/2-1).*drdx;

% d(w*PL) = d(PL)*w + PL*d(w)
pl = spdiags(weight,0,N,N)*dpl + spdiags(dweight,0,N,N)*pl;
