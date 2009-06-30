% Given the degrees of freedom and the domain scaling parameter, returns the
% resulting affine scaling parameter so that all the p/gq nodes lie inside [0,L]
function[scale] = derive_scale(N,L,s,alpha)

pl_parameters;

[x,w] = gq(N,s,alpha);

scale = L/max(x);
