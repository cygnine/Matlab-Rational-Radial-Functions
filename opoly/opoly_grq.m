function[x,w] = opoly_grq(as,bs,N,a);

% function[x,w] = opoly_grq(as,bs,N);
% Returns the N-point Gauss-Radau polynomial quadrature for the orthogonal
% polynomials corresponding to the recurrence coefficients as and bs.
% The fixed Radau point is at the location a
%
% 20080522: acn

as = as(1:N);
bs = bs(1:N);

temp = eval_opoly(a,as,bs,[N-2,N-1]);
as(N) = a - bs(N)*temp(1)/temp(2);

as = as(:);
bs = bs(:);

J = spdiags([sqrt([bs(2:end);0]) as sqrt(bs)], -1:1, N, N);

[v,d] = eig(full(J));
x = diag(d);
w = bs(1)*v(1,:).^2.';
