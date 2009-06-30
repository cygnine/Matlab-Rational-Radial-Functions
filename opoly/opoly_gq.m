function[x,w] = opoly_gq(as,bs,N);

% function[x,w] = opoly_gq(as,bs,N);
% Returns the N-point Gaussian polynomial quadrature for the orthogonal
% polynomials corresponding to the recurrence coefficients as and bs
%
% 20080522: acn

as = as(1:N);
bs = bs(1:N);

as = as(:);
bs = bs(:);

J = spdiags([sqrt([bs(2:end);0]) as sqrt(bs)], -1:1, N, N);

[v,d] = eig(full(J));
x = diag(d);
w = bs(1)*v(1,:).^2.';
